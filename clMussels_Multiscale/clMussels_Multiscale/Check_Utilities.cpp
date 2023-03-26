//
//  Check_Utilities.cpp
//  clMussels_IBM
//
//  Created by Johan Van de Koppel on 29/11/2014.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <math.h>
#include <assert.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif


void CheckSortDevice(cl_command_queue command_queue,
                     cl_mem d_Hash, unsigned int* h_Hash, unsigned int* h_OldHash,
                     cl_mem d_X, float* h_X, float* h_OldX,
                     cl_mem d_Index, unsigned int numElements, int Counter)
{
    
    unsigned int Position_Vector_Memory = sizeof(float) * numElements;
    unsigned int Hash_Vector_Memory = sizeof(unsigned int) * numElements;
    
    // Declaring the host arrays
    unsigned int* h_Index = (unsigned int *)malloc(Hash_Vector_Memory);
    
    // Copy the old h_X and h_Hash values to the backup memory
    for(int i=0;i<numElements;i++){
        h_OldHash[i]=h_Hash[i];
        h_OldX[i]=h_X[i];
    }
    
    size_t ret; // declaring the variable that receives error messages
    
    // Get the new values from dX, d_Hash, and d_Index from the device
    ret = clEnqueueReadBuffer(command_queue, d_X, CL_TRUE, 0, (size_t)(size_t)Position_Vector_Memory, h_X, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, d_Hash, CL_TRUE, 0, (size_t)Hash_Vector_Memory, h_Hash, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, d_Index, CL_TRUE, 0, (size_t)Hash_Vector_Memory, h_Index, 0, NULL, NULL);
    
    
    printf("\n-------------------------------------------\n");
    
    printf("\nStarting iteration %d\n", Counter);
    
    printf("\n h_OldX ; h_OldHash - i  ;  h_X    ; h_Hash - h_Index\n");

    for(int i=0;i<20;i++){
        printf(" %4.3f ; %5d - %5d  ;  %4.3f ; %5d - %5d\n",
               h_OldX[i], h_OldHash[i],i, h_X[i], h_Hash[i],h_Index[i]);
    }
    
    printf("\n Check the order: ");
    
    using namespace std;
    
    // first see if the final list is ordered
    for(unsigned int i=0;i<(numElements-1);i++){
        if (!(h_Hash[i] <= h_Hash[i+1])) {
            cout << " Error at " << i<< " " << h_Hash[i]<<" ,"<<i+1<<" "<<h_Hash[i+1]<<endl;
            exit(10);
        }
    }
    cout << "OK !" << endl;
    
    cout << " Check the permutation: ";
    
    // check if the permutation corresponds to the original list
    for(unsigned int i=0;i<numElements;i++){
        if (!(h_Hash[i] == h_OldHash[h_Index[i]])) {
            cout << " Error at agent "<< i << "!" <<endl;
            exit(10);
        }
    }
    cout << "OK !"<<endl;
    
    cout << " Check the permutation on X: ";
    
    // check if the permutation corresponds to the original list
    for(unsigned int i=0;i<numElements;i++){
        if (!(h_X[i] == h_OldX[h_Index[i]])) {
            cout << " Error at agent "<< i << "!" <<endl;
            exit(10);
        }
    }
    cout << "OK !"<<endl;
}

void CheckHashDevice(cl_command_queue command_queue,
                     cl_mem d_HashBinStart,
                     cl_mem d_HashBinEnd,
                     size_t HashSize)
{
    
    unsigned int* h_HashBinStart = (unsigned int *)malloc((HashSize)*sizeof(unsigned int));
    unsigned int* h_HashBinEnd = (unsigned int *)malloc((HashSize)*sizeof(unsigned int));
    unsigned int Hash_Memory = ((unsigned int)HashSize)*sizeof(unsigned int);
    
    size_t ret;
    
    ret = clEnqueueReadBuffer(command_queue, d_HashBinStart, CL_TRUE, 0, (size_t)Hash_Memory, h_HashBinStart, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, d_HashBinEnd, CL_TRUE, 0, (size_t)Hash_Memory, h_HashBinEnd, 0, NULL, NULL);

    // Check the Hashmapping results
    for(int i=0;i<(HashSize);i++){
        printf(" %d - %d\n", h_HashBinStart[i],h_HashBinEnd[i]);
    }; printf("\n");
}

void CheckRandomDevice(cl_command_queue command_queue,
                     cl_mem R_Angle,
                     cl_mem R_Step,
                     unsigned int numElements)
{
    unsigned int Position_Vector_Memory = sizeof(float) * numElements;
    
    float* h_RandAngle = (float *)malloc((Position_Vector_Memory));
    float* h_RandStep = (float *)malloc((Position_Vector_Memory));
    
    size_t ret;

    // Check the Random Numbers
    ret = clEnqueueReadBuffer(command_queue, R_Angle, CL_TRUE, 0, (size_t)Position_Vector_Memory,
                              h_RandAngle, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, R_Step, CL_TRUE, 0, (size_t)Position_Vector_Memory,
                              h_RandStep, 0, NULL, NULL);

    printf("=> Random Numbers: \n");
    for(int i=0;i<100;i++){
        printf("%4.3f ", h_RandStep[i]);
    }; printf(" \n\n");
}

