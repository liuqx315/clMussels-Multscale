//
//  Mussel_IBM_Kernel.cl
//  TestSort
//
//  Created by Johan Van de Koppel on 21/11/2014.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#include "Settings_and_Parameters.h"

////////////////////////////////////////////////////////////////////////////////
//----- Kernel that reorders the X and Y arrays on basis of sorting Index ------
////////////////////////////////////////////////////////////////////////////////
__kernel void ReorderXY(__global unsigned int* Index,
                        __global float* X,    __global float* Y,
                        __global float* Xnew, __global float* Ynew)
{
    size_t current = get_global_id(0);
    
    Xnew[current] = X[Index[current]];
    Ynew[current] = Y[Index[current]];
}
////////////////////////////////////////////////////////////////////////////////
//----- Kernel that resets the HashBin Start and End values --------------------
////////////////////////////////////////////////////////////////////////////////
__kernel void HashTableReset(__global unsigned int* Start, __global unsigned int* End)
{
    size_t current = get_global_id(0);
    
    Start[current]=-1;
    End[current]=-1;
}

////////////////////////////////////////////////////////////////////////////////
//---------- Finding the first and last mussel index for each hash bin ---------
////////////////////////////////////////////////////////////////////////////////

__kernel void Find_Start_End_Each_Hash(__global unsigned int* d_Hash, __global unsigned int* Start, __global unsigned int* End,
                               const unsigned int IndNr, const unsigned int HS )
{
    size_t current = get_global_id(0);
    int CurrentHash, PreviousHash, NextHash;
    
    CurrentHash=d_Hash[current];
    
    // Makes a table for each hashbin, which are the fist and last mussel
    // Looks for the positions where the hashnr changes, with exceptions for
    // the first and the last.
    if ( (current>0) && (current<(IndNr-1)) ){

        PreviousHash=d_Hash[current-1];
        NextHash=d_Hash[current+1];
        
        if(CurrentHash!=PreviousHash) Start[CurrentHash]=current;
        if(CurrentHash!=NextHash) End[CurrentHash]=current;
    } else {
        if (current==0) // For the first mussel
        {
            NextHash=d_Hash[current+1];
            
            Start[CurrentHash]=0;
            if(CurrentHash!=NextHash) End[CurrentHash]=current;
        }
        if (current==(IndNr-1)) // For the last mussel
        {
            PreviousHash=d_Hash[current-1];
            
            End[CurrentHash]=IndNr-1;
            if(CurrentHash!=PreviousHash) Start[CurrentHash]=current;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//---------- The kernel calculating the movement -------------------------------
////////////////////////////////////////////////////////////////////////////////

__kernel void Calc_Movement (
                    __global float* d_X,     __global float*d_Y,      // The position vectors
                    __global float* d_Xnew,  __global float*d_Ynew,   // The new position vectors
                    __global unsigned int*  d_Index,                          // Mussel index nr
                    __global float* R_Angle, __global float* R_Step,  // Random Numbers
                    __global unsigned int*  d_Hash,                           // The hashmap vectors
                    __global unsigned int*  Begin,   __global unsigned int* End,      // The hashbin vectors
                    const unsigned int HD,                                    // Hash Dimension
                    const unsigned int HS,                                    // Hash Size
                    unsigned int IndNum)                                      // Number of individuals
{
    
    #define D1_2 (D1*D1)
    #define D2_2 (D2*D2)
    
    float diffx,diffy;
    float Distance2, Dens1, Dens2, Beta, Stepsize, Angle;
    int N1=0;
    int N2=0;
    
    size_t current = get_global_id(0);
    
    int CurrentHash=d_Hash[current];
    int CheckHash;
    
    if (current<IndNum)
    {
        // For each hashbin directly around the central center, the distances
        // to other mussels are calculated
        for(int k=-1; k<2; k++) {
            for(int l=-1; l<2; l++) {
                
                CheckHash = CurrentHash + k*HD + l;
                
                // Check if Hash Bin is within the domain and
                // whether it has mussels (if not, bin start = -1)
                if( (CheckHash>=0) && (CheckHash<HS) && (Begin[CheckHash]!=-1)){
                    
                    for(int j=Begin[CheckHash]; j<=End[CheckHash]; j++) {
                        diffx=d_X[j]-d_X[current];
                        diffy=d_Y[j]-d_Y[current];
                        Distance2 = diffx*diffx + diffy*diffy;
                        
                        // If a mussels is closer than D1 or D2, than N1 or N2 are
                        // incremented, respectively
                        N1 = N1 + ( (Distance2 <= (float)D1_2) ? 1 : 0);
                        N2 = N2 + ( (Distance2 <= (float)D2_2) ? 1 : 0);
                    }
                }
            }
        }
        
        Dens1 = (N1-1)/(D1*D1*pi); // Density is mussel number devided by surface
        Dens2 = (N2-1)/(D2*D2*pi); // and again for N2
        
        Angle=(float)R_Angle[current]*2; // pi is omited as we use sinpi & cospi
        
        // Calculating the parameter Beta of the gamma (Exp) distribution
        Beta = 1/(max(0.001,P1*Dens1+P2*Dens2)+P3);
        Stepsize = - Beta * log((float)R_Step[current]);
        
        // Updating the position into a new variable
        d_Xnew[current] = d_X[current] + Stepsize*sinpi(Angle);
        d_Ynew[current] = d_Y[current] + Stepsize*cospi(Angle);
        
        float D_Size = (float)DomainSize;
        
        // Periodic boundary conditions
        if (d_Xnew[current]>D_Size)  { d_Xnew[current] = d_Xnew[current] - D_Size; }
        if (d_Xnew[current]<0.0) { d_Xnew[current] = d_Xnew[current] + D_Size; }
        if (d_Ynew[current]>D_Size)  { d_Ynew[current] = d_Ynew[current] - D_Size; }
        if (d_Ynew[current]<0.0) { d_Ynew[current] = d_Ynew[current] + D_Size; }
        
        // Calculating Hash Dimension (cm per hashbin)
        unsigned int HashX=(unsigned int)(d_X[current]/(DomainSize/HD));
        unsigned int HashY=(unsigned int)(d_Y[current]/(DomainSize/HD));
        d_Hash[current]=HashY*HD+HashX;
        d_Index[current]=current;
    }
    
};



