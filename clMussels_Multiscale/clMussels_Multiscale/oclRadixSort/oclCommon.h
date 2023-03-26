
#ifndef OCLCOMMON_H
#define OCLCOMMON_H

#define check(err) if (err != CL_SUCCESS) printf("line %d, error: %d\n", __LINE__, err); 
#define debug(str) printf("%s\n", str);

#include <sys/time.h>
#include <iostream>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

//timing function. Returns time between consecutive function calls in seconds.
extern double getTime();
extern cl_program getProgramFromKernelFile(std::string filename, cl_context context, cl_int *err);
extern void checkProgramBuild(cl_int errorCode, cl_device_id deviceId, cl_program program);
extern unsigned int roundPow2(unsigned int x);
unsigned int iSnapUp(unsigned int dividend, unsigned int divisor);


#endif
