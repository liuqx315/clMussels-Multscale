/*
* Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
*
* Please refer to the NVIDIA end user license agreement (EULA) associated
* with this source code for terms and conditions that govern your use of
* this software. Any use, reproduction, disclosure, or distribution of
* this software and related documentation outside the terms of the EULA
* is strictly prohibited.
*
*/

//----------------------------------------------------------------------------
// Scans each warp in parallel ("warp-scan"), one element per thread.
// uses 2 numElements of shared memory per thread (64 = elements per warp)
//----------------------------------------------------------------------------
#define WARP_SIZE 32
uint scanwarp(uint val, volatile __local uint* sData, int maxlevel)
{
    // The following is the same as 2 * RadixSort::WARP_SIZE * warpId + threadInWarp = 
    // 64*(threadIdx.x >> 5) + (threadIdx.x & (RadixSort::WARP_SIZE - 1))
    int localId = get_local_id(0);
    int idx = 2 * localId - (localId & (WARP_SIZE - 1));
    sData[idx] = 0;
    idx += WARP_SIZE;
    sData[idx] = val;     

    if (0 <= maxlevel) { sData[idx] += sData[idx - 1]; }
    if (1 <= maxlevel) { sData[idx] += sData[idx - 2]; }
    if (2 <= maxlevel) { sData[idx] += sData[idx - 4]; }
    if (3 <= maxlevel) { sData[idx] += sData[idx - 8]; }
    if (4 <= maxlevel) { sData[idx] += sData[idx -16]; }

    return sData[idx] - val;  // convert inclusive -> exclusive
}

//----------------------------------------------------------------------------
// scan4 scans 4*RadixSort::CTA_SIZE numElements in a block (4 per thread), using 
// a warp-scan algorithm
//----------------------------------------------------------------------------
uint4 scan4(uint4 idata, __local uint* ptr)
{    
    
    uint idx = get_local_id(0);

    uint4 val4 = idata;
    uint sum[3];
    sum[0] = val4.x;
    sum[1] = val4.y + sum[0];
    sum[2] = val4.z + sum[1];
    
    uint val = val4.w + sum[2];
    
    val = scanwarp(val, ptr, 4);
    barrier(CLK_LOCAL_MEM_FENCE);

    if ((idx & (WARP_SIZE - 1)) == WARP_SIZE - 1)
    {
        ptr[idx >> 5] = val + val4.w + sum[2];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

	if (idx < WARP_SIZE)
		ptr[idx] = scanwarp(ptr[idx], ptr, 2);
    
    barrier(CLK_LOCAL_MEM_FENCE);

    val += ptr[idx >> 5];

    val4.x = val;
    val4.y = val + sum[0];
    val4.z = val + sum[1];
    val4.w = val + sum[2];

    return val4;
}

#ifdef MAC
__kernel uint4 rank4(uint4 preds, __local uint* sMem, __local uint* numtrue)
#else
uint4 rank4(uint4 preds, __local uint* sMem, __local uint* numtrue)
#endif
{
	int localId = get_local_id(0);
	int localSize = get_local_size(0);
	uint4 address = scan4(preds, sMem);
	
	if (localId == localSize - 1) {
		numtrue[0] = address.w + preds.w;
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	
	uint4 rank;
	int idx = localId*4;
	rank.x = (preds.x) ? address.x : numtrue[0] + idx - address.x;
	rank.y = (preds.y) ? address.y : numtrue[0] + idx + 1 - address.y;
	rank.z = (preds.z) ? address.z : numtrue[0] + idx + 2 - address.z;
	rank.w = (preds.w) ? address.w : numtrue[0] + idx + 3 - address.w;	
	return rank;
}
__kernel void findRadixOffsets(__global uint2* keys,
				__global uint* counters,
				__global uint* blockOffsets,
				uint startbit,
				uint numElements,
				uint totalBlocks,
				__local uint* sRadix1){
	
	__local uint  sStartPointers[16];

	uint groupId = get_group_id(0);
	uint localId = get_local_id(0);
   	uint groupSize = get_local_size(0);
    	uint2 radix2 = keys[get_global_id(0)];

        
	sRadix1[2 * localId]     = (radix2.x >> startbit) & 0xF;
	sRadix1[2 * localId + 1] = (radix2.y >> startbit) & 0xF;

	// Finds the position where the sRadix1 entries differ and stores start 
	// index for each radix.
	if(localId < 16) {
        	sStartPointers[localId] = 0; 
    	}
    	
	barrier(CLK_LOCAL_MEM_FENCE);

    	if((localId > 0) && (sRadix1[localId] != sRadix1[localId - 1]) ) {
        	sStartPointers[sRadix1[localId]] = localId;
    	}
    	if(sRadix1[localId + groupSize] != sRadix1[localId + groupSize - 1]) {
        	sStartPointers[sRadix1[localId + groupSize]] = localId + groupSize;
    	}
    	
	barrier(CLK_LOCAL_MEM_FENCE);

	if(localId < 16) {
        	blockOffsets[groupId*16 + localId] = sStartPointers[localId];
    	}
    	
	barrier(CLK_LOCAL_MEM_FENCE);

    	// Compute the sizes of each block.
    	if((localId > 0) && (sRadix1[localId] != sRadix1[localId - 1]) ) {
        	sStartPointers[sRadix1[localId - 1]] = localId - sStartPointers[sRadix1[localId - 1]];
    	}
    
	if(sRadix1[localId + groupSize] != sRadix1[localId + groupSize - 1] ) {
        sStartPointers[sRadix1[localId + groupSize - 1]] = localId + groupSize - sStartPointers[sRadix1[localId + groupSize - 1]];}
        
	if(localId == groupSize - 1) {
        	sStartPointers[sRadix1[2 * groupSize - 1]] = 2 * groupSize - sStartPointers[sRadix1[2 * groupSize - 1]];
    	}
    	
	barrier(CLK_LOCAL_MEM_FENCE);

	if(localId < 16) {
        	counters[localId * totalBlocks + groupId] = sStartPointers[localId];
    	}
}

void localSort(uint4 *key, 
			uint4 *values, 
			uint nbits, 
			uint startbit, 
			__local uint* sMem, 
			__local uint* numtrue){
	
	int localId = get_local_id(0);
	int localSize = get_local_size(0);
	
	for(uint shift = startbit; shift < (startbit + nbits); ++shift){
		uint4 lsb;
		lsb.x = !(((*key).x >> shift) & 0x1);
		lsb.y = !(((*key).y >> shift) & 0x1);
        	lsb.z = !(((*key).z >> shift) & 0x1);
        	lsb.w = !(((*key).w >> shift) & 0x1);
        
		uint4 r;
		
		r = rank4(lsb, sMem, numtrue);

        	// This arithmetic strides the ranks across 4 CTA_SIZE regions
        	sMem[(r.x & 3) * localSize + (r.x >> 2)] = (*key).x;
		sMem[(r.y & 3) * localSize + (r.y >> 2)] = (*key).y;
		sMem[(r.z & 3) * localSize + (r.z >> 2)] = (*key).z;
		sMem[(r.w & 3) * localSize + (r.w >> 2)] = (*key).w;
		barrier(CLK_LOCAL_MEM_FENCE);

		// The above allows us to read without 4-way bank conflicts:
		(*key).x = sMem[localId];
		(*key).y = sMem[localId +     localSize];
		(*key).z = sMem[localId + 2 * localSize];
		(*key).w = sMem[localId + 3 * localSize];

		barrier(CLK_LOCAL_MEM_FENCE);

        	// This arithmetic strides the ranks across 4 CTA_SIZE regions
		sMem[(r.x & 3) * localSize + (r.x >> 2)] = (*values).x;
		sMem[(r.y & 3) * localSize + (r.y >> 2)] = (*values).y;
		sMem[(r.z & 3) * localSize + (r.z >> 2)] = (*values).z;
		sMem[(r.w & 3) * localSize + (r.w >> 2)] = (*values).w;
		barrier(CLK_LOCAL_MEM_FENCE);

		// The above allows us to read without 4-way bank conflicts:
		(*values).x = sMem[localId];
		(*values).y = sMem[localId +     localSize];
		(*values).z = sMem[localId + 2 * localSize];
		(*values).w = sMem[localId + 3 * localSize];

		barrier(CLK_LOCAL_MEM_FENCE);
	}
}

__kernel void sortWorkGroup(__global uint4* keysIn, 
				__global uint4* keysOut,
				__global uint4* valuesIn,
				__global uint4* valuesOut,					
				uint nbits,
				uint startbit,
				uint numElements, 
				uint totalBlocks,
				__local uint* sMem){
	
	__local uint numtrue[1];

	int globalId = get_global_id(0);
	uint4 key, value;
	key = keysIn[globalId];
	value = valuesIn[globalId];
	barrier(CLK_LOCAL_MEM_FENCE);
	
	localSort(&key, &value, nbits, startbit, sMem, numtrue);
	
	keysOut[globalId] = key;
	valuesOut[globalId] = value;

}

__kernel void reorderWorkGroup(__global uint  *outKeys, 
                                __global uint2  *keys, 
                                __global uint *outValues,
				__global uint2 *values,
				__global uint  *blockOffsets, 
                                __global uint  *offsets, 
                                __global uint  *sizes, 
                                uint startbit,
                                uint numElements,
                                uint totalBlocks,
                                __local uint2* sKeys2,
				__local uint2* sValues2){
    

	__local uint *sKeys1 = (__local uint*)sKeys2; 
	__local uint *sValues1 = (__local uint*)sValues2; 

	__local uint sOffsets[16];
	__local uint sBlockOffsets[16];

	uint groupId = get_group_id(0);
	uint globalId = get_global_id(0);
	uint localId = get_local_id(0);
	uint groupSize = get_local_size(0);

    	sKeys2[localId]   = keys[globalId];
	sValues2[localId] = values[globalId];

	if(localId < 16){
        	sOffsets[localId]      = offsets[localId * totalBlocks + groupId];
        	sBlockOffsets[localId] = blockOffsets[groupId * 16 + localId];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	uint radix = (sKeys1[localId] >> startbit) & 0xF;
	uint globalOffset = sOffsets[radix] + localId - sBlockOffsets[radix];

	if (globalOffset < numElements){
        	outKeys[globalOffset]   = sKeys1[localId];
        	outValues[globalOffset]   = sValues1[localId];
	}

	radix = (sKeys1[localId + groupSize] >> startbit) & 0xF;
	globalOffset = sOffsets[radix] + localId + groupSize - sBlockOffsets[radix];

	if (globalOffset < numElements){
        	outKeys[globalOffset]   = sKeys1[localId + groupSize];
    		outValues[globalOffset] = sValues1[localId+groupSize];
	}
}

