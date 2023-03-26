//
//  IBM_program.h
//  clMussels_Multiscale
//
//  Created by Johan Van de Koppel on 23/12/2014.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#ifndef IBM_PROGRAM_H
#define IBM_PROGRAM_H

#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <vector>

using namespace std;

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "Settings_and_Parameters.h"
#include "Forward_Definitions.h"

#include "oclRadixSort.h"
#include "oclCommon.h"

class IBM_Mussels{
    
private:
    
    unsigned int numIndividuals;
    unsigned int numElements;
    static const unsigned int maxElements = Maximum_Nr_Elements;
    static const unsigned int minElements = Minimum_Nr_Elements;
    
    cl_mem d_X;
    cl_mem d_Y;
    
    cl_mem d_Xnew;
    cl_mem d_Ynew;
    
    cl_mem d_Index;
    cl_mem d_Hash;
    
    cl_mem d_HashBinStart;
    cl_mem d_HashBinEnd;
    
    // The array that receives the random numbers
    cl_mem R_Angle;
    cl_mem R_Step;
    
    cl_mem d_MT;
    unsigned int nPerRng;

    cl_command_queue IBM_command_queue;
    cl_context IBM_context;
    cl_program IBM_program;
    cl_program MT_program;
    cl_device_id* devices;
    
    cl_kernel ReorderXYkernel;
    cl_kernel HashTableResetKernel;
    cl_kernel Find_Start_End_kernel;
    cl_kernel MovementKernel;
    
    cl_kernel MersenneTwisterKernel;
    oclRadixSort* fSort;
    
    size_t global_item_size;
    size_t local_item_size;
    
    size_t global_Hash_size;
    size_t local_Hash_size;
    
    unsigned int HashDimension;
    unsigned int HashSize;
    
    unsigned int Position_Vector_Memory;
    unsigned int Position_Used_Memory;
    unsigned int Hash_Vector_Memory;
    unsigned int Hash_Used_Memory;
    unsigned int Hash_Memory;
    
    void ReorderXY();
    void MakeHashTable();
    void Find_Start_End();
    void Movement();
    void RandomNumers_MT();
    void UpdateCoordinates();
    
    ProfileClassCPU* GPUProfiler;
    
public:
    IBM_Mussels(cl_device_id deviceId, cl_context context, cl_command_queue queue);
    ~IBM_Mussels();
    
    int GetHashDimension();
    
    void RunSimulation(vector<Individual> & Mussel, unsigned int numIndividuals,
                       unsigned int MovementCycles, ProfileClassCPU & Profiler);
};

#endif

