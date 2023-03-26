
#include "oclScan.h"
#include "oclCommon.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

oclScan::oclScan(cl_device_id deviceId, cl_context context, cl_command_queue queue, unsigned int numElements){

	
	cl_int err;

	fNumElements = numElements;

	if (numElements > fMaxInclusiveScanSize){
		fTempBuffer = clCreateBuffer(context, 
					CL_MEM_READ_WRITE, 
					numElements / fMaxInclusiveScanSize * sizeof(cl_uint), 
					NULL, 
					&err);
		check(err);
	}

	//assign pointers
	fContext = context;
	fQueue = queue;

	fProgram = getProgramFromKernelFile("oclScan.cl", fContext, &err);
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

	//Create kernel objects
	fScanExclusiveLocal1Kernel = clCreateKernel(fProgram, "scanExclusiveLocal1", &err);
	fScanExclusiveLocal2Kernel = clCreateKernel(fProgram, "scanExclusiveLocal2", &err);
	fScanExclusiveSmallKernel = clCreateKernel(fProgram, "scanExclusiveSmall", &err);

	fUniformUpdateKernel = clCreateKernel(fProgram, "uniformUpdate", &err);
	fNaiveKernel = clCreateKernel(fProgram, "scanNaive", &err);

	size_t Multiple;
	err =  clGetKernelWorkGroupInfo(fNaiveKernel, deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &Multiple, 	NULL);
	


	check(err);
}

oclScan::~oclScan(){
	cl_int err;

	err  = clReleaseKernel(fScanExclusiveLocal1Kernel);
	err |= clReleaseKernel(fScanExclusiveLocal2Kernel);
	err |= clReleaseKernel(fUniformUpdateKernel);
	err |= clReleaseKernel(fNaiveKernel);
	err |= clReleaseKernel(fScanExclusiveSmallKernel);


	if (fNumElements > fMaxInclusiveScanSize)
    		err |= clReleaseMemObject(fTempBuffer);
    	err |= clReleaseProgram(fProgram);
	check(err);
}

void oclScan::scan(cl_mem d_Dst, cl_mem d_Src, unsigned int batchSize, unsigned int arrayLength){
		
	//printf("arraylength is %u\n", arrayLength);

	if ((arrayLength >= 4) && (arrayLength <= 256)){
		scanNaive(d_Dst, d_Src, arrayLength);	
		//printf("naive scan\n");
	} else if ( (arrayLength>256) &&(arrayLength <= 1024)){
		//need this to sort bins 16384 
		scanExclusiveShort(d_Dst, d_Src, batchSize, arrayLength);	
		//printf("exclusive short\n");
	} else if (arrayLength > 1024) {
		scanExclusiveLarge(d_Dst, d_Src, batchSize, arrayLength);
		//printf("exclusive large\n");	
	}
}

void oclScan::scanNaive(cl_mem d_Dst, cl_mem d_Src, unsigned int arrayLength){

	size_t g[1] = {arrayLength};
	size_t l[1] = {arrayLength};
	cl_uint err;	
	unsigned int extraSpace = (arrayLength / 16);
	unsigned int sharedMem = sizeof(unsigned int)*extraSpace;

	err  = clSetKernelArg(fNaiveKernel, 0, sizeof(cl_mem), (void *)&d_Dst);
	err |= clSetKernelArg(fNaiveKernel, 1, sizeof(cl_mem), (void *)&d_Src);
	err |= clSetKernelArg(fNaiveKernel, 2, sizeof(unsigned int), (void *)&arrayLength);
	err |= clSetKernelArg(fNaiveKernel, 3, 2*sharedMem, NULL);

	check(err);
	err = clEnqueueNDRangeKernel(fQueue, fNaiveKernel, 1, NULL, g, l, 0, NULL, NULL);
    	check(err);

}

void oclScan::scanExclusiveShort(cl_mem d_Dst, cl_mem d_Src, unsigned int batchSize, unsigned int arrayLength){
  
	cl_int err;
    	size_t l[1] = {fLocalSize};
	size_t g[1] = {arrayLength/2};

	
	err  = clSetKernelArg(fScanExclusiveLocal1Kernel, 0, sizeof(cl_mem), (void *)&d_Dst);
	err |= clSetKernelArg(fScanExclusiveLocal1Kernel, 1, sizeof(cl_mem), (void *)&d_Src);
	err |= clSetKernelArg(fScanExclusiveLocal1Kernel, 2, 2*fLocalSize * sizeof(unsigned int), NULL);
	err |= clSetKernelArg(fScanExclusiveLocal1Kernel, 3, sizeof(unsigned int), (void *)&arrayLength);
    	check(err);


	err = clEnqueueNDRangeKernel(fQueue, fScanExclusiveLocal1Kernel, 1, NULL, g, l, 0, NULL, NULL);
    	check(err);

    return;
}

// main exclusive scan routine
void oclScan::scanExclusiveLarge(cl_mem d_Dst, cl_mem d_Src, unsigned int batchSize, unsigned int arrayLength){
    

	scanExclusiveLocal1(d_Dst,
        			d_Src,
        			(batchSize * arrayLength) / (4 * fLocalSize),
        			4 * fLocalSize);

    	scanExclusiveLocal2(fTempBuffer,
        			d_Dst,
        			d_Src,
        			batchSize,
        			arrayLength / (4 * fLocalSize));

	uniformUpdate(d_Dst,
        		fTempBuffer,
        		(batchSize * arrayLength) / (4 * fLocalSize));

}

void oclScan::scanExclusiveLocal1(cl_mem d_Dst,
    					cl_mem d_Src,
    					unsigned int n,
    					unsigned int size){
   	
	cl_int err;
    	size_t l[1] = {fLocalSize};
	size_t g[1] = {(n*size)/4};

	err  = clSetKernelArg(fScanExclusiveLocal1Kernel, 0, sizeof(cl_mem), (void *)&d_Dst);
	err |= clSetKernelArg(fScanExclusiveLocal1Kernel, 1, sizeof(cl_mem), (void *)&d_Src);
	err |= clSetKernelArg(fScanExclusiveLocal1Kernel, 2, 2 * fLocalSize * sizeof(unsigned int), NULL);
	err |= clSetKernelArg(fScanExclusiveLocal1Kernel, 3, sizeof(unsigned int), (void *)&size);
    	check(err);


	err = clEnqueueNDRangeKernel(fQueue, fScanExclusiveLocal1Kernel, 1, NULL, g, l, 0, NULL, NULL);
    	check(err);

}

void oclScan::scanExclusiveLocal2(cl_mem d_Buffer,
    				cl_mem d_Dst,
				cl_mem d_Src,
    				unsigned int n,
    				unsigned int size){
	cl_int err;
	unsigned int elements = n * size;
	size_t l[1] = {fLocalSize};
	size_t g[1] = {iSnapUp(elements, fLocalSize)};

	
	err  = clSetKernelArg(fScanExclusiveLocal2Kernel, 0, sizeof(cl_mem), (void *)&d_Buffer);
	err |= clSetKernelArg(fScanExclusiveLocal2Kernel, 1, sizeof(cl_mem), (void *)&d_Dst);
	err |= clSetKernelArg(fScanExclusiveLocal2Kernel, 2, sizeof(cl_mem), (void *)&d_Src);
	err |= clSetKernelArg(fScanExclusiveLocal2Kernel, 3, 2 * fLocalSize * sizeof(unsigned int), NULL);
	

	err |= clSetKernelArg(fScanExclusiveLocal2Kernel, 4, sizeof(unsigned int), (void *)&elements);
	err |= clSetKernelArg(fScanExclusiveLocal2Kernel, 5, sizeof(unsigned int), (void *)&size);
	check(err);

	
     	err = clEnqueueNDRangeKernel(fQueue, fScanExclusiveLocal2Kernel, 1, NULL, g, l, 0, NULL, NULL);
    	check(err);


}

void oclScan::uniformUpdate(cl_mem d_Dst,
    			cl_mem d_Buffer,
    			unsigned int n){
	cl_int err;
     	size_t l[1] = {fLocalSize};
    	size_t g[1] = {fLocalSize*n};



	err  = clSetKernelArg(fUniformUpdateKernel, 0, sizeof(cl_mem), (void *)&d_Dst);
	err |= clSetKernelArg(fUniformUpdateKernel, 1, sizeof(cl_mem), (void *)&d_Buffer);
	check(err);
    
	err = clEnqueueNDRangeKernel(fQueue, fUniformUpdateKernel, 1, NULL, g, l, 0, NULL, NULL);
    	check(err);

}
