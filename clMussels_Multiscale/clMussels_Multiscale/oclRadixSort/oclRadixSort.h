#ifndef OCLRADIXSORT_H
#define OCLRADIXSORT_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "oclScan.h"

class oclRadixSort{

private:

	//Radix Sort Buffers
	cl_mem fTempKeysBuffer;
	cl_mem fTempValuesBuffer;
	cl_mem fCountersBuffer;
	cl_mem fCountersSumBuffer;
	cl_mem fOffsetsBuffer;
	
	cl_command_queue fQueue;
	cl_context fContext;
	cl_program fProgram;	

	//Radix Sort Kernels
	cl_kernel fSortWorkGroupKernel;
	cl_kernel fFindRadixOffsetsKernel;
	cl_kernel fReorderWorkGroupKernel;
	
	oclScan *fScan;

	static const unsigned int fLocalSize = 256; 
	static const unsigned int fBitStep = 4;
	static const unsigned int fWarpSize = 32;

	void sortWorkGroup(cl_mem d_keys, cl_mem d_values,  unsigned int nbits, unsigned int startbit, unsigned int numElements);
	void getOffsets(unsigned int startbit, unsigned int numElements);
	void reorderWorkGroup(cl_mem d_keys, cl_mem d_values, unsigned int startbit, unsigned int numElements);

public:
	oclRadixSort(cl_device_id device_id, cl_context context, cl_command_queue queue, unsigned int numElements);
	~oclRadixSort();

	void sort(cl_mem d_keys, cl_mem d_values, unsigned int numElements, unsigned int keyBits);
	void step(cl_mem keys, cl_mem values, unsigned int startbit, unsigned int nbits);

	void sortSingleBlock();

};

#endif
