//
//  Forward_Definitions.h
//  clMussels_Multiscale
//
//  Created by Johan Van de Koppel on 22/12/2014.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#ifndef FORWARD_DEFINITIONS_H
#define FORWARD_DEFINITIONS_H

#include "Profiler.h"

struct Individual  // Structure defining a mussel
{
    float x;
    float y;
    int GridNr;
};

struct HashBinStruct   // Structure defining a HashBin
{
    int Start;
    int End;
};

////////////////////////////////////////////////////////////////////////////////
// Forward definitions from functions at the end of this code file
////////////////////////////////////////////////////////////////////////////////

// Forward definitions for "PDE_Computing_Kernel.cpp"
void MusselsKernel ( float*, float*, float*, int );

// Forward definitions for "Transfer_between_PDE_IBM.cpp"
void AdaptMusselNr(float*, vector<Individual> &, int &, ProfileClassCPU &);
void AdaptMusselGrid(float* &, vector<Individual> , int, ProfileClassCPU &);

// Forward definitions for "HostFunctions.cpp"
void Save_PDE_Matrices(char*, float*, float*, float*, int);
void Save_XY_Vectors(char*, std::vector<Individual> &, int, int);

// Forward definitions for "Initializations.cpp"
void randomInit (float*, int, int, int);
void InitMusselsNr(vector<Individual> &, float *, int &);

#endif
