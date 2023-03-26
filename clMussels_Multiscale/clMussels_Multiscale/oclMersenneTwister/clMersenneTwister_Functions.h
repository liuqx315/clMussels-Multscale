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

#include <stdio.h>
#include <math.h>
#include <iostream>

#include "clMersenneTwister.h"

//#define UINT32_C(value) \
//   __CONCAT(value, UL)

static mt_struct MT[MT_RNG_COUNT];
//static uint32_t state[MT_NN];

///////////////////////////////////////////////////////////////////////////////
//Load twister configurations
///////////////////////////////////////////////////////////////////////////////
void loadMTGPU(const char *fname,
               const unsigned int seed,
               mt_struct_stripped *h_MT,
               const size_t size)
{
    FILE* fd = 0;
    #ifdef _WIN32
    // open the file for binary read
        errno_t err;
        if ((err = fopen_s(&fd, fname, "rb")) != 0)
    #else
        // open the file for binary read
        if ((fd = fopen(fname, "rb")) == 0)
    #endif
        {
            if(fd)
            {
                fclose (fd);
            }
        }
    
    for (unsigned int i = 0; i < size; i++)
        fread(&h_MT[i], sizeof(mt_struct_stripped), 1, fd);
    fclose(fd);
    
    for(unsigned int i = 0; i < size; i++)
        h_MT[i].seed = rand(); // JvdK: changed from "seed" to avoid banding
}

extern "C" void initMTRef(const char *fname){
    
    FILE* fd = 0;
    #ifdef _WIN32
        // open the file for binary read
        errno_t err;
        if ((err = fopen_s(&fd, fname, "rb")) != 0)
    #else
        // open the file for binary read
        if ((fd = fopen(fname, "rb")) == 0)
    #endif
        {
            if(fd)
            {
                fclose (fd);
            }
        }

    for (int i = 0; i < MT_RNG_COUNT; i++)
    {
        //Inline structure size for compatibility,
        //since pointer types are 8-byte on 64-bit systems (unused *state variable)
        fread(MT + i, 16 /* sizeof(mt_struct) */ * sizeof(int), 1, fd);
    }

    fclose(fd);
}



void BoxMuller(float& u1, float& u2) {
    float   r = sqrtf(-2.0f * logf(u1));
    float phi = 2 * PI * u2;
    u1 = r * cosf(phi);
    u2 = r * sinf(phi);
}

//Initialize/seed twister for current GPU context
void RandomInfo(size_t N_per_Range, int Gr_Height, int Gr_Width){
	printf(" A Mersenne Twister - based Random number generator\n");
	printf(" for use on Graphics Processors. Based on the code \n");
	printf(" by Makoto Matsumoto and Takuji Nishimura (2001). \n\n");
	
	printf(" Dimensions   : %d x %d\n", Gr_Height,Gr_Width);
	printf(" Grid Size    : %d\n", Gr_Height*Gr_Width);
 	printf(" MT_RNG_COUNT : %d\n", MT_RNG_COUNT);
	printf(" N_PER_RNG    : %d\n", (int) N_per_Range);
	printf(" RAND_N       : %d\n\n", MT_RNG_COUNT * (int)N_per_Range);
    
}

/*
 * Victor Podlozhnyuk @ 05/13/2007:
 * 1) Fixed sgenrand_mt():
 *    - Fixed loop indexing: 'i' variable was off by one.
 *    - apply wmask right on the state element initialization instead
 *      of separate loop, which could produce machine-dependent results(wrong).
 * 2) Slightly reformatted sources to be included into CUDA SDK.
 */

#define SHIFT1 18

void sgenrand_mt(uint32_t seed, mt_struct *mts){
    int i;
    
    mts->state[0] = seed & mts->wmask;
    
    for(i = 1; i < mts->nn; i++){
        mts->state[i] = (UINT32_C(1812433253) * (mts->state[i - 1] ^ (mts->state[i - 1] >> 30)) + i) & mts->wmask;
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
    }
    mts->i = mts->nn;
}

uint32_t genrand_mt(mt_struct *mts){
    uint32_t *st, uuu, lll, aa, x;
    int k,n,m,lim;
    
    if(mts->i >= mts->nn ){
        n = mts->nn; m = mts->mm;
        aa = mts->aaa;
        st = mts->state;
        uuu = mts->umask; lll = mts->lmask;
        
        lim = n - m;
        for(k = 0; k < lim; k++){
            x = (st[k]&uuu)|(st[k+1]&lll);
            st[k] = st[k + m] ^ (x >> 1) ^ (x&1U ? aa : 0U);
        }
        
        lim = n - 1;
        for(; k < lim; k++){
            x = (st[k] & uuu)|(st[k + 1] & lll);
            st[k] = st[k + m - n] ^ (x >> 1) ^ (x & 1U ? aa : 0U);
        }
        
        x = (st[n - 1] & uuu)|(st[0] & lll);
        st[n - 1] = st[m - 1] ^ (x >> 1) ^ (x&1U ? aa : 0U);
        mts->i=0;
    }
    
    x = mts->state[mts->i];
    mts->i += 1;
    x ^= x >> mts->shift0;
    x ^= (x << mts->shiftB) & mts->maskB;
    x ^= (x << mts->shiftC) & mts->maskC;
    x ^= x >> mts->shift1;
    
    return x;
}


