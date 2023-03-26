//
//  IBM_program.cpp
//  clMussels_Multiscale
//
//  Created by Johan Van de Koppel on 23/12/2014.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#include "IBM_program.h"

#include "clMersenneTwister.h"
#include "clMersenneTwister_Functions.h"

#include "Settings_and_Parameters.h"

cl_program BuildKernelFile(std::string, cl_context, cl_device_id*, cl_int*);

// Forward definitions from functions at the end of this code file
void CheckSort(unsigned int* , unsigned int* , float*, float*, unsigned int*, unsigned int);
void CheckSortDevice(cl_command_queue,cl_mem,unsigned int*,unsigned int*,cl_mem,float*,float*,cl_mem,unsigned int, int);
void CheckHashDevice(cl_command_queue,cl_mem,cl_mem,size_t);
void CheckRandomDevice(cl_command_queue,cl_mem,cl_mem,unsigned int);

IBM_Mussels::IBM_Mussels(cl_device_id deviceId,
                           cl_context context,
                           cl_command_queue queue)
{
    //Constructor, build stuff in here
    cl_int err;
    IBM_context = context;
    IBM_command_queue = queue;
    numElements = 0;
    numIndividuals = 0;
    
    //----------Constant and variable definition--------------------------------
    
    // The location of the code is obtain from the __FILE__ macro
    const std::string SourcePath (__FILE__);
    const std::string PathName = SourcePath.substr (0,SourcePath.find_last_of("/")+1);
    
    const std::string DataName = "Mussels.dat";
    const std::string PRNGDataName = "MT_Output.dat";
    
    const std::string IncludeFolder = "-I " + PathName;
    const std::string DataPath = PathName + DataName;
    const std::string PRNGDataPath = PathName + DataName;
    
    // ---- Calculating hash dimensions ----------------------------------------
    // Calculating the HashBinSize as the first value above D2 that fits within
    // the domain without creating a fraction.
    
    unsigned int HBS = D2;   // HBS = Hash Bin Size
    while ((DomainSize%HBS)!=0) { HBS++; };
    HashDimension = (unsigned int)DomainSize/HBS;
    HashSize = HashDimension*HashDimension;
    
    Position_Vector_Memory = sizeof(float) * maxElements;
    Hash_Vector_Memory = sizeof(unsigned int) * maxElements;
    Hash_Memory = ((unsigned int)HashSize)*sizeof(unsigned int);
    
    //----------Linking the Kernel files----------------------------------------
    IBM_program = BuildKernelFile("IBM_Computing_Kernel.cl", context,
                                            &deviceId, &err);
    
    MT_program = BuildKernelFile("oclMersenneTwister/clMersenneTwister.cl", context,
                                            &deviceId, &err);
    
    //----------Create memory buffers-------------------------------------------
    d_X = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t)Position_Vector_Memory, NULL, &err);
    d_Y = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t)Position_Vector_Memory, NULL, &err);
    
    d_Xnew = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t)Position_Vector_Memory, NULL, &err);
    d_Ynew = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t)Position_Vector_Memory, NULL, &err);
    
    d_Index = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t)Hash_Vector_Memory, NULL, &err);
    d_Hash  = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t)Hash_Vector_Memory, NULL, &err);
    
    d_HashBinStart  = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t)Hash_Memory, NULL, &err);
    d_HashBinEnd  = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t)Hash_Memory, NULL, &err);
    
    // The array that receives the random numbers
    R_Angle = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t) Position_Vector_Memory, NULL, &err);
    R_Step = clCreateBuffer(IBM_context, CL_MEM_READ_WRITE, (size_t) Position_Vector_Memory, NULL, &err);

    check(err);
    
    //----------Create ReorderXY kernel-----------------------------------------
    ReorderXYkernel = clCreateKernel(IBM_program, "ReorderXY", &err);
    if (err!=0) printf(" > Create ReorderXY Kernel Error number: %d \n\n", err);

    //----------Reset the Hash Tabel--------------------------------------------
    HashTableResetKernel = clCreateKernel(IBM_program, "HashTableReset", &err);
    if (err!=0) printf(" > Create HashTable Reset Kernel Error number: %d \n\n", err);

    //----------Create Find_Start_End_Each_Hash kernel--------------------------
    Find_Start_End_kernel = clCreateKernel(IBM_program, "Find_Start_End_Each_Hash", &err);
    if (err!=0) printf(" > Create Find_Start_End_kernel Error number: %d \n\n", err);
    
    //----------Create Movement Kernel------------------------------------------
    MovementKernel = clCreateKernel(IBM_program, "Calc_Movement", &err);
    if (err!=0) printf(" > Create MovementKernel Error number: %d \n\n", err);
    
    //----------Mersenne Twister specific settings and variables----------------
    
    cl_int ciErr1, ciErr2;                          // Error code var
    const int seed = 66e6;                          // Seeding number for Random Number generation (not used)
    nPerRng = maxElements/MT_RNG_COUNT;   // # of recurrence steps, must be even if do Box-Muller transformation
    //const int nRand = MT_RNG_COUNT * nPerRng;       // Output size
    
    const std::string MTDataName = PathName + "oclMersenneTwister/MT_data/MersenneTwister.dat";
    const std::string MTRawName = PathName + "oclMersenneTwister/MT_data/MersenneTwister.raw";
    const std::string MTIncludeFolder = "-I " + PathName + "oclMersenneTwister/";
    
    mt_struct_stripped* h_MT = (mt_struct_stripped*)malloc(sizeof(mt_struct_stripped)*MT_RNG_COUNT);
    MersenneTwisterKernel = NULL;             // Mersenne Twister kernel
    
    loadMTGPU(MTDataName.c_str(), seed, h_MT, MT_RNG_COUNT);
    initMTRef(MTRawName.c_str());
    
    // Setting up the MT algorithm data on the GPU
    d_MT = clCreateBuffer(IBM_context, CL_MEM_READ_ONLY, sizeof(mt_struct_stripped)*MT_RNG_COUNT, NULL, &ciErr2);
    if(ciErr2!=0) {printf("Error allocating memory to the device");}
    ciErr1 |= clEnqueueWriteBuffer(IBM_command_queue, d_MT, CL_TRUE, 0,
                                   sizeof(mt_struct_stripped)*MT_RNG_COUNT, h_MT, 0, NULL, NULL);
    
    //----------Linking the Mersenne Twister kernel code------------------------
    
    // Linking to the specific kernel function
    MersenneTwisterKernel = clCreateKernel(MT_program, "MersenneTwister", &ciErr1);
    if (ciErr1!=0) { printf(" > Failed to create the MT kernel: %d \n\n", ciErr1);};

    free(h_MT);

    //----------Set up the sorting process--------------------------------------
    
    fSort = new oclRadixSort(deviceId, IBM_context, IBM_command_queue, maxElements);
    
    //----------Set up the GPU block size and number----------------------------
    
    //global_item_size = numElements;
    //local_item_size = 1024;
    
    //global_Hash_size = HashSize;
    //local_Hash_size = 1;
    
    //----------Set up the Profiling--------------------------------------------
    
    GPUProfiler = new ProfileClassCPU (7, Profiling);
    
}

////////////////////////////////////////////////////////////////////////////////
// The reorder function
////////////////////////////////////////////////////////////////////////////////

void IBM_Mussels::ReorderXY()
{
    int err;
    
    GPUProfiler -> Mark();
    
    /* Set OpenCL kernel arguments */
    err = clSetKernelArg(ReorderXYkernel, 0, sizeof(cl_mem), (void *)&d_Index);
    err = clSetKernelArg(ReorderXYkernel, 1, sizeof(cl_mem), (void *)&d_X);
    err = clSetKernelArg(ReorderXYkernel, 2, sizeof(cl_mem), (void *)&d_Y);
    err = clSetKernelArg(ReorderXYkernel, 3, sizeof(cl_mem), (void *)&d_Xnew);
    err = clSetKernelArg(ReorderXYkernel, 4, sizeof(cl_mem), (void *)&d_Ynew);
    
    /* Execute ReorderXY kernel */
    err = clEnqueueNDRangeKernel(IBM_command_queue, ReorderXYkernel, 1, NULL,
                                 &global_item_size, &local_item_size, 0, NULL, NULL);
    if (err!=0) { printf(" ReorderXYkernel Error number: %d \n\n", err); exit(-10); }
    
    // Copy the Xnew/Ynew to the X and Y vectors
    err = clEnqueueCopyBuffer (IBM_command_queue, d_Xnew, d_X, 0, 0, (size_t)Position_Vector_Memory, 0, NULL, NULL);
    err = clEnqueueCopyBuffer (IBM_command_queue, d_Ynew, d_Y, 0, 0, (size_t)Position_Vector_Memory, 0, NULL, NULL);
    
    if(Profiling==ON) {
        err = clFinish(IBM_command_queue);
        GPUProfiler -> Register(2);
    }

    
}

////////////////////////////////////////////////////////////////////////////////
// The MakeHashTable function
////////////////////////////////////////////////////////////////////////////////

void IBM_Mussels::MakeHashTable()
{
    int err;
 
    GPUProfiler -> Mark();
    
    /* Set OpenCL kernel arguments */
    err  = clSetKernelArg(HashTableResetKernel, 0, sizeof(cl_mem), (void *)&d_HashBinStart);
    err |= clSetKernelArg(HashTableResetKernel, 1, sizeof(cl_mem), (void *)&d_HashBinEnd);
    
    /* Execute HashTableReset */
    err = clEnqueueNDRangeKernel(IBM_command_queue, HashTableResetKernel, 1, NULL,
                                 &global_Hash_size, &local_Hash_size, 0, NULL, NULL);
    if (err!=0) { printf(" > HashTableResetKernel Error number: %d \n\n", err); exit(-10); }
    
    if(Profiling==ON) {
        err = clFinish(IBM_command_queue);
        GPUProfiler -> Register(3);
    }
    
}

////////////////////////////////////////////////////////////////////////////////
// Finding the first and last mussel for each hashbin
////////////////////////////////////////////////////////////////////////////////

void IBM_Mussels::Find_Start_End()
{
    int err;
    
    GPUProfiler -> Mark();
    
    /* Set OpenCL kernel arguments */
    err  = clSetKernelArg(Find_Start_End_kernel, 0, sizeof(cl_mem), (void *)&d_Hash);
    err |= clSetKernelArg(Find_Start_End_kernel, 1, sizeof(cl_mem), (void *)&d_HashBinStart);
    err |= clSetKernelArg(Find_Start_End_kernel, 2, sizeof(cl_mem), (void *)&d_HashBinEnd);
    err |= clSetKernelArg(Find_Start_End_kernel, 3, sizeof(unsigned int), (void *)&numIndividuals);
    err |= clSetKernelArg(Find_Start_End_kernel, 4, sizeof(unsigned int), (void *)&HashSize);
    
    //-------------- Making a Hask Table for distance calculation ------
    err |= clEnqueueNDRangeKernel(IBM_command_queue, Find_Start_End_kernel, 1, NULL,
                                 &global_item_size, &local_item_size, 0, NULL, NULL);
    
    if (err!=0) { printf(" > Find_Start_End_kernel Error number: %d \n\n", err); exit(-10); }
    
    if(Profiling==ON) {
        err = clFinish(IBM_command_queue);
        GPUProfiler -> Register(4);
    }
    
    if(CheckHashing==ON)
        CheckHashDevice(IBM_command_queue,d_HashBinStart,d_HashBinEnd,HashSize);
}

////////////////////////////////////////////////////////////////////////////////
// Makeing a new set of Random numbers for R_Angle and R_Step
////////////////////////////////////////////////////////////////////////////////

void IBM_Mussels::RandomNumers_MT()
{
    GPUProfiler -> Mark();
    
    //-------------- Random Number Generation --------------------------
    
    // Setting up MT Kernel execution parameters
    size_t MTglobalWorkSize= MT_RNG_COUNT;    // 1D var for Total # of work items
    size_t MTlocalWorkSize = nPerRng;         // 1D var for # of work items in the work group
    
    cl_int err;
    
    // Setting the parameters of the kernel function
    err  = clSetKernelArg(MersenneTwisterKernel, 0, sizeof(cl_mem), (void*)&R_Angle);
    err |= clSetKernelArg(MersenneTwisterKernel, 1, sizeof(cl_mem), (void*)&d_MT);
    err |= clSetKernelArg(MersenneTwisterKernel, 2, sizeof(int),    (void*)&nPerRng);
    
    err |= clEnqueueNDRangeKernel(IBM_command_queue, MersenneTwisterKernel, 1, NULL,
                                     &MTglobalWorkSize, &MTlocalWorkSize, 0, NULL, NULL);
    
    err |= clSetKernelArg(MersenneTwisterKernel, 0, sizeof(cl_mem), (void*)&R_Step);
    
    err |= clEnqueueNDRangeKernel(IBM_command_queue, MersenneTwisterKernel, 1, NULL,
                                     &MTglobalWorkSize, &MTlocalWorkSize, 0, NULL, NULL);
    
    if(Profiling==ON) {
        err = clFinish(IBM_command_queue);
        GPUProfiler -> Register(5);
    }
    
    check(err);

}

////////////////////////////////////////////////////////////////////////////////
// Calculating the movement of each mussel
////////////////////////////////////////////////////////////////////////////////

void IBM_Mussels::Movement()
{
    GPUProfiler -> Mark();
    
    //-------------- Calculating movement ------------------------------
    
    cl_int err;
    
    /* Set OpenCL kernel arguments */
    err  = clSetKernelArg(MovementKernel,  0, sizeof(cl_mem), (void *)&d_X);
    err |= clSetKernelArg(MovementKernel,  1, sizeof(cl_mem), (void *)&d_Y);
    err |= clSetKernelArg(MovementKernel,  2, sizeof(cl_mem), (void *)&d_Xnew);
    err |= clSetKernelArg(MovementKernel,  3, sizeof(cl_mem), (void *)&d_Ynew);
    err |= clSetKernelArg(MovementKernel,  4, sizeof(cl_mem), (void *)&d_Index);
    err |= clSetKernelArg(MovementKernel,  5, sizeof(cl_mem), (void *)&R_Angle);
    err |= clSetKernelArg(MovementKernel,  6, sizeof(cl_mem), (void *)&R_Step);
    err |= clSetKernelArg(MovementKernel,  7, sizeof(cl_mem), (void *)&d_Hash);
    err |= clSetKernelArg(MovementKernel,  8, sizeof(cl_mem), (void *)&d_HashBinStart);
    err |= clSetKernelArg(MovementKernel,  9, sizeof(cl_mem), (void *)&d_HashBinEnd);
    err |= clSetKernelArg(MovementKernel, 10, sizeof(unsigned int),   (void *)&HashDimension);
    err |= clSetKernelArg(MovementKernel, 11, sizeof(unsigned int),   (void *)&HashSize);
    err |= clSetKernelArg(MovementKernel, 12, sizeof(unsigned int),   (void *)&numIndividuals);
    
    err |= clEnqueueNDRangeKernel(IBM_command_queue, MovementKernel, 1, NULL,
                                 &global_item_size, &local_item_size, 0, NULL, NULL);
    if (err!=0) { printf(" > Movement Kernel Error number: %d \n\n", err); exit(-10); }
    
    if(Profiling==ON) {
        err = clFinish(IBM_command_queue);
        GPUProfiler -> Register(6);
    }
    
    if(CheckPRNG==ON) CheckRandomDevice(IBM_command_queue,R_Angle,R_Step,numElements);

}

////////////////////////////////////////////////////////////////////////////////
// Update the mussel coordinates from the new values in d_Xnew and d_Ynew
////////////////////////////////////////////////////////////////////////////////

void IBM_Mussels::UpdateCoordinates()
{
    GPUProfiler -> Mark();
    
    //-------------- Update mussel coordinates -------------------------
    
    cl_int err;
    
    // Copy the Xnew/Ynew to the X and Y vectors
    err  = clEnqueueCopyBuffer (IBM_command_queue, d_Xnew, d_X, 0, 0, (size_t)Position_Used_Memory, 0, NULL, NULL);
    err |= clEnqueueCopyBuffer (IBM_command_queue, d_Ynew, d_Y, 0, 0, (size_t)Position_Used_Memory, 0, NULL, NULL);
    
    //err = clFinish(IBM_command_queue);
    
    if(Profiling==ON) {
        err = clFinish(IBM_command_queue);
        GPUProfiler -> Register(7);
    }

    check(err);
}

////////////////////////////////////////////////////////////////////////////////
// Run the Simulation for "MovementCycles" times
////////////////////////////////////////////////////////////////////////////////

void IBM_Mussels::RunSimulation(vector<Individual> & Mussel, unsigned int numInd,
                                unsigned int MovementCycles, ProfileClassCPU & Profiler)
{
    cl_int err;
    
    numIndividuals =numInd;
    numElements = (unsigned int) max(powf(2, (float)trunc(log2((float) numIndividuals)+1)),minElements);
    
    Position_Used_Memory = numElements*sizeof(float);
    Hash_Used_Memory = numElements*sizeof(unsigned int);
    
    // Marking the start for the Profiler, from Main Program
    Profiler.Mark();
    
    if(numElements>maxElements)
    {
        printf(" Program error: mussels number exceeding maximum of %d",
               maxElements);
        exit(10);
    }
    
    //----------Set up the Work group size and number --------------------------
    global_item_size = numElements;
    local_item_size = 512;
    
    global_Hash_size = HashSize;
    local_Hash_size = 1;

    //----------Host variable declaration and allocation -----------------------   
    float* h_X    = (float *)malloc(Position_Used_Memory);
    float* h_OldX = (float *)malloc(Position_Used_Memory);
    float* h_Y    = (float *)malloc(Position_Used_Memory);
    unsigned int*  h_Hash = (unsigned int *) malloc(Hash_Used_Memory);
    unsigned int*  h_OldHash = (unsigned int *)malloc(Hash_Used_Memory);
    unsigned int*  h_Index = (unsigned int *)malloc(Hash_Used_Memory);
    
    // ---- Initialization of the variables ------------------------------------
    
    unsigned int HashX,HashY;
    
    for(int i=0; i<numIndividuals; i++) {
        h_X[i]=Mussel[i].x;
        h_Y[i]=Mussel[i].y;
        h_Index[i]=i;
        
        // Calculating Hash Dimension (cm per hashbin)
        HashX=(unsigned int)h_X[i]/(DomainSize/HashDimension);
        HashY=(unsigned int)h_Y[i]/(DomainSize/HashDimension);
        h_Hash[i]=HashY*HashDimension+HashX;
        h_OldHash[i]=h_Hash[i];
        h_OldX[i]=h_X[i];
    }
    
    for(int i=numIndividuals; i<numElements; i++) {
        h_X[i]=20000000;
        h_Y[i]=20000000;
        h_Index[i]=i;
        
        // Calculating Hash Dimension (cm per hashbin)
        h_Hash[i]=HashSize + 1;
        h_OldHash[i]=h_Hash[i];
        h_OldX[i]=h_X[i];
    }
    
    err = clEnqueueWriteBuffer(IBM_command_queue, d_X, CL_TRUE, 0, (size_t)Position_Used_Memory, h_X, 0, NULL, NULL);
    err = clEnqueueWriteBuffer(IBM_command_queue, d_Y, CL_TRUE, 0, (size_t)Position_Used_Memory, h_Y, 0, NULL, NULL);
    
    err = clEnqueueWriteBuffer(IBM_command_queue, d_Index, CL_TRUE, 0, (size_t)Hash_Used_Memory, h_Index, 0, NULL, NULL);
    err = clEnqueueWriteBuffer(IBM_command_queue, d_Hash, CL_TRUE, 0, (size_t)Hash_Used_Memory, h_Hash, 0, NULL, NULL);

    //-------------- Kernel initiation -----------------------------------------
    
    for (int Counter=0;Counter<MovementCycles;Counter++)
    {
        GPUProfiler -> Mark();
        fSort -> sort(d_Hash,d_Index, numElements, 32);
        if(Profiling==ON) {
            err = clFinish(IBM_command_queue);
            GPUProfiler -> Register(1);
        }
        
        ReorderXY();
        
        if(CheckSorting==ON)
        {  CheckSortDevice(IBM_command_queue, d_Hash, h_Hash, h_OldHash,
                           d_X, h_X, h_OldX, d_Index, numIndividuals, Counter); }
        
        MakeHashTable();
        Find_Start_End();
        RandomNumers_MT();
        Movement();
        UpdateCoordinates();
        
        if(CheckSorting==ON)
        {
            err  = clEnqueueReadBuffer(IBM_command_queue, d_X, CL_TRUE, 0, (size_t)Position_Used_Memory, h_X, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(IBM_command_queue, d_Y, CL_TRUE, 0, (size_t)Position_Used_Memory, h_Y, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(IBM_command_queue, d_Hash, CL_TRUE, 0, (size_t)Hash_Used_Memory, h_Hash, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(IBM_command_queue, d_Index, CL_TRUE, 0, (size_t)Hash_Used_Memory, h_Index, 0, NULL, NULL);
            
            check(err);
        }
    }

    err  = clEnqueueReadBuffer(IBM_command_queue, d_X, CL_TRUE, 0, (size_t)Position_Used_Memory, h_X, 0, NULL, NULL);
    err |= clEnqueueReadBuffer(IBM_command_queue, d_Y, CL_TRUE, 0, (size_t)Position_Used_Memory, h_Y, 0, NULL, NULL);
    err |= clEnqueueReadBuffer(IBM_command_queue, d_Hash, CL_TRUE, 0, (size_t)Hash_Used_Memory, h_Hash, 0, NULL, NULL);
    
    check(err);
    
    err = clFlush(IBM_command_queue);
    err = clFinish(IBM_command_queue);
    
    for(int i=0; i<numIndividuals; i++) {
        Mussel[i].x=h_X[i];
        Mussel[i].y=h_Y[i];
    }
    
    // Freeing host space
    free(h_X);
    free(h_OldX);
    free(h_Y);
    
    free(h_Hash);
    free(h_OldHash);
    free(h_Index);
    
    // Registering the time spend
    Profiler.Register(3);

}

////////////////////////////////////////////////////////////////////////////////
//----------------------Hash Dimension Reporting function-----------------------
////////////////////////////////////////////////////////////////////////////////

int IBM_Mussels::GetHashDimension(){
    return HashDimension;
}

////////////////////////////////////////////////////////////////////////////////
//----------------------Clean up memory-----------------------------------------
////////////////////////////////////////////////////////////////////////////////


IBM_Mussels::~IBM_Mussels()
{
    
    //----------Report on GPU profiling-----------------------------------------
    if(Profiling==ON) {
        printf(" GPU TIME USAGE \n\n");
        printf(" Sorting took                  : %4.3f seconds\n", GPUProfiler -> Report(1));
        printf(" Reordering took               : %4.3f seconds\n", GPUProfiler -> Report(2));
        printf(" Cleaning the Hashtable took   : %4.3f seconds\n", GPUProfiler -> Report(3));
        printf(" Finding start and end took    : %4.3f seconds\n", GPUProfiler -> Report(4));
        printf(" Random number generation took : %4.3f seconds\n", GPUProfiler -> Report(5));
        printf(" Calculating movement took     : %4.3f seconds\n", GPUProfiler -> Report(6));
        printf(" Updating positions took       : %4.3f seconds\n", GPUProfiler -> Report(7));
    }
    
    cl_int err;
    
    // Freeing kernel and block space
    err = clFlush(IBM_command_queue);
    err = clFinish(IBM_command_queue);
    
    // Freeing device variables
    err = clReleaseMemObject(d_X);
    err = clReleaseMemObject(d_Y);
    err = clReleaseMemObject(d_Xnew);
    err = clReleaseMemObject(d_Ynew);
    
    err = clReleaseMemObject(d_Index);
    err = clReleaseMemObject(d_Hash);
    
    err = clReleaseMemObject(d_HashBinStart);
    err = clReleaseMemObject(d_HashBinEnd);
    
    err = clReleaseMemObject(R_Step);
    err = clReleaseMemObject(R_Angle);
    
    err = clReleaseKernel(ReorderXYkernel);
    err = clReleaseKernel(HashTableResetKernel);
    err = clReleaseKernel(Find_Start_End_kernel);
    err = clReleaseKernel(MovementKernel);
    err = clReleaseKernel(MersenneTwisterKernel);
    
    err = clReleaseProgram(IBM_program);
    err = clReleaseProgram(MT_program);
    err = clReleaseMemObject(d_MT);
    
    delete(fSort);
    delete(GPUProfiler);

}


