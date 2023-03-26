#include <iostream>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "oclRadixSort.h"
#include "oclCommon.h"

oclRadixSort::oclRadixSort(cl_device_id deviceId, 
				cl_context context, 
				cl_command_queue queue, 
				unsigned int numElements){
	
	//Constructor, build stuff in here
	cl_int err;
	fContext = context;
	fQueue = queue;

	unsigned int numWorkGroups = ((numElements % (4*fLocalSize)) == 0) ? 
            (numElements / (4*fLocalSize)) : (numElements / (4*fLocalSize) + 1);

	fProgram = getProgramFromKernelFile("oclRadixSort.cl", fContext, &err);
	check(err);

	//Build Program
	
    err = clBuildProgram(fProgram, 0, NULL, "-cl-fast-relaxed-math", NULL, NULL);
    if (err != CL_SUCCESS){
    //return details of any compilation error
        size_t len;
        char buffer[2048];
        clGetProgramBuildInfo(fProgram, deviceId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("error is %d, buffer is %s\n", err, buffer);
        return;
   	}

	//create memory buffers
	fTempKeysBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(unsigned int) * numElements, NULL, &err);
	fTempValuesBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * numElements, NULL, &err);
	fCountersBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, fWarpSize * numWorkGroups * sizeof(unsigned int), NULL, &err);
	fCountersSumBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, fWarpSize * numWorkGroups * sizeof(unsigned int), NULL, &err);
	fOffsetsBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, fWarpSize * numWorkGroups * sizeof(unsigned int), NULL, &err); 

	check(err);


	//build kernels
	fSortWorkGroupKernel = clCreateKernel(fProgram, "sortWorkGroup", &err);
	fFindRadixOffsetsKernel = clCreateKernel(fProgram, "findRadixOffsets", &err);
	fReorderWorkGroupKernel = clCreateKernel(fProgram, "reorderWorkGroup", &err);
	check(err);

	//create instance of oclScan
	fScan = new oclScan(deviceId, context, queue, numElements/2/fLocalSize*16);


}

oclRadixSort::~oclRadixSort(){
	//free memory buffers	
	clReleaseMemObject(fTempKeysBuffer);
	clReleaseMemObject(fTempValuesBuffer);
	clReleaseMemObject(fCountersBuffer);
	clReleaseMemObject(fCountersSumBuffer);
	clReleaseMemObject(fOffsetsBuffer);

	clReleaseProgram(fProgram);

	//free kernels
	clReleaseKernel(fSortWorkGroupKernel);
	clReleaseKernel(fFindRadixOffsetsKernel);
	clReleaseKernel(fReorderWorkGroupKernel);
	
	delete fScan; 

}

void oclRadixSort::sort(cl_mem keys, cl_mem values, unsigned int numElements, unsigned int keyBits){
	
	int i = 0;
	unsigned int len = numElements/2/fLocalSize*16;
	
	while (keyBits > i*fBitStep){

		sortWorkGroup(keys, values, fBitStep, i*fBitStep, numElements);
		getOffsets(i*fBitStep, numElements);	
		fScan->scan(fCountersSumBuffer, fCountersBuffer, 1, len);
		reorderWorkGroup(keys, values, i*fBitStep, numElements);
		i++;
	}

}

void oclRadixSort::sortWorkGroup(cl_mem keys, cl_mem values,  unsigned int nbits, unsigned int startbit, unsigned int numElements){
	

	unsigned int totalWorkGroups = numElements/4/fLocalSize;
	size_t g[1] = {fLocalSize*totalWorkGroups};
	size_t l[1] = {fLocalSize};
	cl_int err;

	//set kernel arguments
	err  = clSetKernelArg(fSortWorkGroupKernel, 0, sizeof(cl_mem), (void*)&keys);
    err |= clSetKernelArg(fSortWorkGroupKernel, 1, sizeof(cl_mem), (void*)&fTempKeysBuffer);
	err |= clSetKernelArg(fSortWorkGroupKernel, 2, sizeof(cl_mem), (void*)&values);
    err |= clSetKernelArg(fSortWorkGroupKernel, 3, sizeof(cl_mem), (void*)&fTempValuesBuffer);
	err |= clSetKernelArg(fSortWorkGroupKernel, 4, sizeof(unsigned int), (void*)&nbits);
	err |= clSetKernelArg(fSortWorkGroupKernel, 5, sizeof(unsigned int), (void*)&startbit);
    err |= clSetKernelArg(fSortWorkGroupKernel, 6, sizeof(unsigned int), (void*)&numElements);
    err |= clSetKernelArg(fSortWorkGroupKernel, 7, sizeof(unsigned int), (void*)&totalWorkGroups);
    err |= clSetKernelArg(fSortWorkGroupKernel, 8, 4*fLocalSize*sizeof(unsigned int), NULL);
    
	//execute kernel 
    err |= clEnqueueNDRangeKernel(fQueue, fSortWorkGroupKernel, 1, NULL, g, l, 0, NULL, NULL);
	check(err);
}

void oclRadixSort::getOffsets(unsigned int startbit, unsigned int numElements){

	unsigned int numWorkGroups = numElements/2/fLocalSize;
	size_t g[1] = {fLocalSize*numWorkGroups};
	size_t l[1] = {fLocalSize};
	cl_int err;
	
	err  = clSetKernelArg(fFindRadixOffsetsKernel, 0, sizeof(cl_mem), (void*)&fTempKeysBuffer);
	err |= clSetKernelArg(fFindRadixOffsetsKernel, 1, sizeof(cl_mem), (void*)&fCountersBuffer);
    err |= clSetKernelArg(fFindRadixOffsetsKernel, 2, sizeof(cl_mem), (void*)&fOffsetsBuffer);
	err |= clSetKernelArg(fFindRadixOffsetsKernel, 3, sizeof(unsigned int), (void*)&startbit);
	err |= clSetKernelArg(fFindRadixOffsetsKernel, 4, sizeof(unsigned int), (void*)&numElements);
	err |= clSetKernelArg(fFindRadixOffsetsKernel, 5, sizeof(unsigned int), (void*)&numWorkGroups);
	err |= clSetKernelArg(fFindRadixOffsetsKernel, 6, 2 * fLocalSize *sizeof(unsigned int), NULL);
	err |= clEnqueueNDRangeKernel(fQueue, fFindRadixOffsetsKernel, 1, NULL, g, l, 0, NULL, NULL);
	check(err);
}
	
void oclRadixSort::reorderWorkGroup(cl_mem keys, cl_mem values, unsigned int startbit, unsigned int numElements){
	
	unsigned int numWorkGroups = numElements/2/fLocalSize;
	size_t g[1] = {fLocalSize*numWorkGroups};
	size_t l[1] = {fLocalSize};
	cl_int err;
	err  = clSetKernelArg(fReorderWorkGroupKernel, 0, sizeof(cl_mem), (void*)&keys);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 1, sizeof(cl_mem), (void*)&fTempKeysBuffer);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 2, sizeof(cl_mem), (void*)&values);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 3, sizeof(cl_mem), (void*)&fTempValuesBuffer);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 4, sizeof(cl_mem), (void*)&fOffsetsBuffer);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 5, sizeof(cl_mem), (void*)&fCountersSumBuffer);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 6, sizeof(cl_mem), (void*)&fCountersBuffer);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 7, sizeof(unsigned int), (void*)&startbit);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 8, sizeof(unsigned int), (void*)&numElements);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 9, sizeof(unsigned int), (void*)&numWorkGroups);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 10, 2 * fLocalSize * sizeof(unsigned int), NULL);
	err |= clSetKernelArg(fReorderWorkGroupKernel, 11, 2 * fLocalSize * sizeof(unsigned int), NULL);	
	err |= clEnqueueNDRangeKernel(fQueue, fReorderWorkGroupKernel, 1, NULL, g, l, 0, NULL, NULL);
	check(err);

}
	
