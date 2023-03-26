//
//  Main_Program.cpp
//
//  Created by Johan Van de Koppel on 03-09-14.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <math.h>
#include <vector>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

using namespace std;

#include "Settings_and_Parameters.h"
#include "Device_Utilities.h"
#include "IBM_program.h"

////////////////////////////////////////////////////////////////////////////////
// Main program code for Mussels
////////////////////////////////////////////////////////////////////////////////

int main()
{
    
    /*----------Defining and allocating memeory on host-----------------------*/
    
    // Defining and allocating the memory blocks for P, W, and O on the host (h)
    float* h_A = (float *)malloc(Grid_Width*Grid_Height*sizeof(float));
	float* h_M = (float *)malloc(Grid_Width*Grid_Height*sizeof(float));
	float* h_S = (float *)malloc(Grid_Width*Grid_Height*sizeof(float));
    
    /*----------Initializing the host arrays----------------------------------*/
    
    srand(40); // Seeding the random number generator
    
	randomInit(h_A, Grid_Width, Grid_Height, ALGAE);
	randomInit(h_M, Grid_Width, Grid_Height, MUSSELS);
    randomInit(h_S, Grid_Width, Grid_Height, SEDIMENT);
    
    /*----------Printing info to the screen ----------------------------------*/

	//system("clear");
    printf("\n");
	printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n");
	printf(" * Mussel bed Patterns at multiple spatial scale         * \n");
	printf(" * OpenCL implementation : Johan van de Koppel, 2014     * \n");
	printf(" * Following a model by Liu et al 2014                   * \n");
	printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n");
    
    /*----------Setting up the devices and loading the kernels ---------------*/
    
    cl_device_id* devices;
    cl_int err;
    
    cl_context context = CreateGPUcontext(devices);
    
    // Print the name of the device that is used
    printf(" Implementing PDE on device %d: ", PDE_Device_No);
    print_device_info(devices, (int)PDE_Device_No);
           
    printf(" Implementing IBM on device %d: ", IBM_Device_No);
    print_device_info(devices, (int)IBM_Device_No);
    
    printf("\n");
    
    // Create the PDE and IBM command queues on the devices
    cl_command_queue PDE_command_queue = clCreateCommandQueue(context, devices[PDE_Device_No], 0, &err);
    cl_command_queue IBM_command_queue = clCreateCommandQueue(context, devices[IBM_Device_No], 0, &err);
    
    /*----------Defining and allocating PDE memeory on device-----------------*/
    
    cl_mem d_A = clCreateBuffer(context, CL_MEM_READ_WRITE, Grid_Memory, NULL, &err);
    cl_mem d_M = clCreateBuffer(context, CL_MEM_READ_WRITE, Grid_Memory, NULL, &err);
    cl_mem d_S = clCreateBuffer(context, CL_MEM_READ_WRITE, Grid_Memory, NULL, &err);
    
    /* Copy input data to the memory buffer */
    err = clEnqueueWriteBuffer(PDE_command_queue, d_A, CL_TRUE, 0, Grid_Width*Grid_Height*sizeof(float), h_A, 0, NULL, NULL);
    err = clEnqueueWriteBuffer(PDE_command_queue, d_M, CL_TRUE, 0, Grid_Width*Grid_Height*sizeof(float), h_M, 0, NULL, NULL);
    err = clEnqueueWriteBuffer(PDE_command_queue, d_S, CL_TRUE, 0, Grid_Width*Grid_Height*sizeof(float), h_S, 0, NULL, NULL);

    /*----------Building the PDE kernel---------------------------------------*/
    
    cl_program PDE_program = BuildKernelFile("PDE_Computing_Kernel.cl", context, &devices[PDE_Device_No], &err);
    if (err!=0)  printf(" > Compile Program Error number: %d \n\n", err);
    
    /* Create OpenCL kernel */
    cl_kernel PDE_kernel = clCreateKernel(PDE_program, "MusselsKernel", &err);
    if (err!=0) printf(" > Create Kernel Error number: %d \n\n", err);
    
    /* Set OpenCL kernel arguments */
    err = clSetKernelArg(PDE_kernel, 0, sizeof(cl_mem), (void *)&d_A);
    err = clSetKernelArg(PDE_kernel, 1, sizeof(cl_mem), (void *)&d_M);
    err = clSetKernelArg(PDE_kernel, 2, sizeof(cl_mem), (void *)&d_S);
    
    /*----------Defining the individual mussels on the host-------------------*/
     
     // Make a host vector that holds the mussels
     vector<Individual> Mussel;
     int MusselNr;
     InitMusselsNr(Mussel, h_M, MusselNr);  //not really needed.
    
    /*----------Building the IBM kernels--------------------------------------*/
    
    // Allocates an object of class "IBM_Mussels", defined in "IBM_program.h"
    IBM_Mussels* IBM = new IBM_Mussels(devices[IBM_Device_No], context, IBM_command_queue);
    
    /*----------Preparing the simulation--------------------------------------*/
    
    float EndTime_EX =(float)EndTime/(float)EX;
    int MovementCycles=(int)(EndTime_EX/(float)NumFrames*(float)CyclesPerTimestep);
    
    // fprintf(stderr," MovementCycles : %d\n\n", MovementCycles);
    
    size_t PDE_global_item_size[] = {Grid_Width*Grid_Height};
    size_t PDE_local_item_size[] = {Block_Size_X*Block_Size_Y};
    
    // -------------- Profiling initialization ---------------------------------
    
    ProfileClassCPU Profiler (4, Profiling);
    
    //----------Report on setup-------------------------------------------------
    
    printf(" Current grid dimensions: %d x %d cells\n",
           Grid_Width, Grid_Height);
    
    printf(" Domain length: %d cm, with %dx%d  hash bins\n",
           DomainSize, IBM->GetHashDimension(),IBM->GetHashDimension());
    
    printf(" Maximal mussel number: %d , initial number: %d\n\n",
           Maximum_Nr_Elements, MusselNr);
    
    /* Progress bar initiation */
    int RealBarWidth=min((int)NumFrames,(int)ProgressBarWidth);
    int BarCounter=0;
    float BarThresholds[RealBarWidth];
    for (int i=0;i<RealBarWidth;i++) {BarThresholds[i] = (float)(i+1)/RealBarWidth*NumFrames;};
    
    if(CheckNumbers!=ON)
    {
        /* Print the reference bar */
        printf(" Progress: [");
        for (int i=0;i<RealBarWidth;i++) { printf("-"); }
        printf("]\n");
        fprintf(stderr, "           >");
    }
    
    // -------------- Kernel initiation ----------------------------------------

	for (int Counter=0;Counter<NumFrames;Counter++)
    {
        // Setting a starting mark for profiling the PDE
        Profiler.Mark();
        
        for (int Runtime=0;Runtime<(int)(EndTime_EX/(float)NumFrames/dT);Runtime++)
        {
            /* Execute OpenCL kernel as data parallel */
            err= clEnqueueNDRangeKernel(PDE_command_queue, PDE_kernel, 1, NULL,
                                         PDE_global_item_size, PDE_local_item_size, 0, NULL, NULL);
            
            if (err!=0) { printf(" > PDE kernel error number: %d \n\n", err); exit(10);}
        }
        
        /* Transfer result to host */
        err = clEnqueueReadBuffer(PDE_command_queue, d_A, CL_TRUE, 0, Grid_Width*Grid_Height*sizeof(float), h_A, 0, NULL, NULL);
        err = clEnqueueReadBuffer(PDE_command_queue, d_M, CL_TRUE, 0, Grid_Width*Grid_Height*sizeof(float), h_M, 0, NULL, NULL);
        err = clEnqueueReadBuffer(PDE_command_queue, d_S, CL_TRUE, 0, Grid_Width*Grid_Height*sizeof(float), h_S, 0, NULL, NULL);
        
        // Registring the time used for the PDE
        Profiler.Register(1);
        
        AdaptMusselNr(h_M, Mussel, MusselNr, Profiler);
        
        IBM->RunSimulation(Mussel, MusselNr, MovementCycles, Profiler);
        
        AdaptMusselGrid(h_M, Mussel, MusselNr, Profiler);
        
        Save_XY_Vectors((char*)"Data/IBM", Mussel,MusselNr,Counter);
        Save_PDE_Matrices((char*)"Data/PDE",h_A,h_M,h_S,Counter);
        
        if(CheckNumbers==ON)
        {
            // Report the timestep and the number of individuals
            fprintf(stderr,"\n Timestep %1.0f of %1.0f, Nr of individuals : %d",
                   (float)Counter/(float)NumFrames*(float)EndTime/24+1,
                   (float)EndTime/24, MusselNr);
        } else if ((float)(Counter+1)>=BarThresholds[BarCounter])
            {
                // Progress the progress bar if time
                fprintf(stderr,"*");
                BarCounter = BarCounter+1;
            }
    }
    
    fprintf(stderr,"<\n\n");
    
    if(Profiling==ON) {
        printf(" CPU/GPU TIME USAGE \n\n");
        printf(" PDE model (GPU) took    : %4.3f seconds\n", Profiler.Report(1));
        printf(" Adapting mussel nr took : %4.3f seconds\n", Profiler.Report(2));
        printf(" IBM model (GPU) took    : %4.3f seconds\n", Profiler.Report(3));
        printf(" Adapting grid took      : %4.3f seconds\n", Profiler.Report(4));
    }
    
    /*---------------------Report on time spending----------------------------*/

    double Timespanned = Profiler.TimeSpanned();
    
    if(Timespanned<60)
        { printf(" Total processing time   : %4.3f (s) \n\n", Timespanned); }
    else
        { printf(" Total processing time   : %d:%d  \n\n",
                 (int)Timespanned/60, (int)Timespanned%60); }

	/*---------------------Clean up memory------------------------------------*/
	
    // Freeing host space
    free(h_A);
	free(h_M);
	free(h_S);
    
    err = clReleaseMemObject(d_A);
    err = clReleaseMemObject(d_M);
    err = clReleaseMemObject(d_S);
    
    err = clReleaseKernel(PDE_kernel);
    err = clReleaseProgram(PDE_program);
    
    //delete(IBM);
    
    err = clReleaseCommandQueue(PDE_command_queue);
    err = clReleaseContext(context);
    
    #if defined(__APPLE__) && defined(__MACH__)
        system("say Simulation finished");
    #endif

	printf("\r Simulation finished! \n\n");
    
}
