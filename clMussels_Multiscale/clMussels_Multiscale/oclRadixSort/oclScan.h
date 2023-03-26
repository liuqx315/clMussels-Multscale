#ifndef OCLSCAN_H
#define OCLSCAN_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

class oclScan {

private:

	unsigned int fNumElements;

	cl_context fContext;
	cl_command_queue fQueue;
	cl_program fProgram;

	//Scan Kernels
	cl_kernel fScanExclusiveLocal1Kernel;
    cl_kernel fScanExclusiveLocal2Kernel;
    cl_kernel fUniformUpdateKernel;
	cl_kernel fNaiveKernel;
	cl_kernel fScanExclusiveSmallKernel;


	cl_mem fTempBuffer;

	static const unsigned int fLocalSize = 256;
	static const unsigned int fMaxInclusiveScanSize = 1024;


	void scanExclusiveLocal1(cl_mem d_Dst, cl_mem d_Src, unsigned int n, unsigned int size);
	void scanExclusiveLocal2(cl_mem d_Buffer, cl_mem d_Dst, cl_mem d_Src, unsigned int n, unsigned int size);
	void uniformUpdate(cl_mem d_Dst, cl_mem d_Buffer, unsigned int n);
	
public:
	oclScan(cl_device_id device_id, cl_context context, cl_command_queue queue, unsigned int numElements);
	~oclScan();
	void scanExclusiveLarge(cl_mem d_Dst, cl_mem d_Src, unsigned int batchSize, unsigned int arrayLength);
	void scan(cl_mem d_Dst, cl_mem d_Src, unsigned int batchSize, unsigned int arrayLength);
	void scanExclusiveShort(cl_mem d_Dst, cl_mem d_Src, unsigned int batchSize, unsigned int arrayLength);
	void scanNaive(cl_mem d_Dst, cl_mem d_Src, unsigned int arrayLength);


};

#endif
