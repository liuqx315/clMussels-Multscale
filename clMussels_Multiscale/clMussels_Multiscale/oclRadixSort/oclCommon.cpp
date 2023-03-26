
#include "oclCommon.h"

#include <stdio.h>
#include <iostream>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

double getTime(){

	static double prev_time = -1;
	struct timeval new_time;
	double result;

	gettimeofday(&new_time, NULL);
	if (prev_time < 0) 
		prev_time = (new_time.tv_usec*1.0e-6)+new_time.tv_sec;

	result = (new_time.tv_usec*1.0e-6)+new_time.tv_sec-prev_time;
	prev_time = (new_time.tv_usec*1.0e-6)+new_time.tv_sec;

	return result;
}

unsigned int iSnapUp(unsigned int dividend, unsigned int divisor){
	return ((dividend % divisor) == 0) ? dividend : (dividend - dividend % divisor + divisor);
}

cl_program getProgramFromKernelFile(std::string filename, cl_context context, cl_int *err){

    // The location of the code is obtain from the __FILE__ macro
    const std::string SourcePath (__FILE__);
    const std::string PathName = SourcePath.substr (0,SourcePath.find_last_of("/")+1);
    
    const std::string FilePath = PathName + filename;
    
    //Read Kernel from file
	FILE *f = fopen(FilePath.c_str(), "r");
	fseek(f, 0, SEEK_END);
	size_t fileSize = ftell(f);
	rewind(f);
	char *fileString = (char*)malloc(fileSize + 1);
	fileString[fileSize] = '\0';
	fread(fileString, sizeof(char), fileSize, f);	
	fclose(f);
	cl_program prog = clCreateProgramWithSource(context, 1, (const char**)&fileString, &fileSize, err);
	free(fileString);
	
	return prog;
} 

void checkProgramBuild(cl_int errorCode, cl_device_id deviceId, cl_program program){

	if (errorCode != CL_SUCCESS){
		//return details of any compilation error
        	size_t len;
        	char buffer[2048];
        	clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        	printf("error is %d, buffer is %s\n", errorCode, buffer);
   	}
}

unsigned int roundPow2(unsigned int x){
	x--;
	x = (x >> 1) | x;
	x = (x >> 2) | x;
	x = (x >> 4) | x;
	x = (x >> 8) | x;
	x = (x >> 16) | x;
	x++; 
	return x;
}




