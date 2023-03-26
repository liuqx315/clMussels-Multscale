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

#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H
#ifndef mersennetwister_h
#define mersennetwister_h

#define ON  1
#define OFF 0

#define AUTORESEED    ON

#define      DCMT_SEED 4172
#define MAX_SOURCE_SIZE (0x100000)

typedef unsigned int uint32_t;
// #define UINT32_C(a) ((uint32_t)a)

typedef struct {
    uint32_t aaa;
    int mm,nn,rr,ww;
    uint32_t wmask,umask,lmask;
    int shift0, shift1, shiftB, shiftC;
    uint32_t maskB, maskC;
    int i;
    uint32_t *state;
}mt_struct;

typedef struct{
  unsigned int matrix_a;
  unsigned int mask_b;
  unsigned int mask_c;
  unsigned int seed;
} mt_struct_stripped;


#define   MT_RNG_COUNT 4096
#define          MT_MM 9
#define          MT_NN 19
#define       MT_WMASK 0xFFFFFFFFU
#define       MT_UMASK 0xFFFFFFFEU
#define       MT_LMASK 0x1U
#define      MT_SHIFT0 12
#define      MT_SHIFTB 7
#define      MT_SHIFTC 15
#define      MT_SHIFT1 18
#define PI 3.14159265358979f

////////////////////////////////////////////////////////////////////////////////
// Prototypes for original Matsumoto & Nishimura code
////////////////////////////////////////////////////////////////////////////////
#if defined(__cplusplus)
extern "C"
#endif
void sgenrand_mt(uint32_t seed, mt_struct *mts);
#if defined(__cplusplus)
extern "C"
#endif
uint32_t genrand_mt(mt_struct *mts);

#endif
#endif
