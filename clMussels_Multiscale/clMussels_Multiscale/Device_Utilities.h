//
//  DeviceUtils.cpp
//
//  Created by Johan Van de Koppel on 03-09-14.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#ifndef DEVICE_UTILITIES_H
#define DEVICE_UTILITIES_H

#include <string>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define ON 1
#define OFF 0

#define Print_All_Devices OFF

void print_platform_info(cl_platform_id platform);
void print_device_info(cl_device_id* device, int Device_Nr);
void Query(cl_uint deviceCount, cl_device_id* devices);
void Get_Build_Errors(cl_program program, cl_device_id* device, cl_int ret);
cl_program BuildKernelFile(std::string filename, cl_context context,
                           cl_device_id* device, cl_int *err);
cl_context CreateGPUcontext(cl_device_id* &devices);

#endif


