//
//  Initializations.cpp
//  cpuMussels_MultiScale
//
//  Created by Johan Van de Koppel on 22/12/2014.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#include <cstdlib>
#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

#include "Settings_and_Parameters.h"
#include "Forward_Definitions.h"

////////////////////////////////////////////////////////////////////////////////
// Allocates a matrix with random float entries
////////////////////////////////////////////////////////////////////////////////

void randomInit (float* data, int x_siz, int y_siz, int type)
{
    int i,j;
    for(i=0;i<y_siz;i++)
    {
        for(j=0;j<x_siz;j++)
        {
            //assigning the first row last row and
            //first column last column as zeroes
            
            if(i==0||i==y_siz-1||j==0||j==x_siz-1)
                data[i*y_siz+j]=0.5f;
            else
            {
                //for every other element find the correct initial
                //value using the conditions below
                if(type==MUSSELS)
                {
                    //printf(" %4.5f ",(rand() / (float)RAND_MAX));
                    if((rand() / (float)RAND_MAX)<0.10f)
                        data[i*y_siz+j] = (float)300.0f;
                    else
                        data[i*y_siz+j] = (float)150.0f;
                    
                }
                else if(type==ALGAE)
                    data[i*y_siz+j]=0.6f;
                else if(type==SEDIMENT)
                    data[i*y_siz+j]=0.0f;
            }
        }
    }
} // End randomInit

////////////////////////////////////////////////////////////////////////////////
// Initialization of the mussel field &Mussel with MusselNo individuals
////////////////////////////////////////////////////////////////////////////////

void InitMusselsNr(vector<Individual> &Mussel, float * h_M, int &MusselNo)
{
    Individual NewMussel;
    int GridX,GridY;
    
    // Checking how many mussels should be in the field
    float M_Sum = 0;
    for(int i=0;i<Grid_Size;i++) M_Sum=M_Sum+h_M[i]*dX*dY;
    int TargetMusselNr = M_Sum*Mussels_Per_Gr;
    
    for(int i=0; i<TargetMusselNr; i++) {
        NewMussel.x=(float)rand()/(float)RAND_MAX*DomainSize;
        NewMussel.y=(float)rand()/(float)RAND_MAX*DomainSize;
        GridX=(int)NewMussel.x/(DomainSize/Grid_Width);
        GridY=(int)NewMussel.y/(DomainSize/Grid_Height);
        NewMussel.GridNr=GridY*Grid_Width+GridX;
        Mussel.push_back(NewMussel);
    }
    MusselNo = (int)Mussel.size();
}

