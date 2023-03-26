//
//  Transfer_between_PDE_IBM.cpp
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
#include <algorithm>

using namespace std;

#include "Settings_and_Parameters.h"
#include "Forward_Definitions.h"

// ---- Sort Container by name function, for the sort function -----------------
bool sortByGrid(const Individual &lhs, const Individual &rhs)
{ return lhs.GridNr < rhs.GridNr; }

////////////////////////////////////////////////////////////////////////////////
// Adjusting the local mussels density to the changes in the grid
////////////////////////////////////////////////////////////////////////////////

void AdaptMusselNr(float* h_M, vector<Individual> &Mussel, int &MusselNo, ProfileClassCPU &Profiler)
{
    int MusselInCell, MusselChange;
    vector<int> MusselAlive (Maximum_Nr_Elements,OFF);
    Individual NewMussel;
    vector<Individual>OldMusselStack,NewMusselStack;
    int GridX,GridY;
    int GridBin[Grid_Height*Grid_Width];
    int CurrentGridNr;
    int counter = 0;
    
    // checking for PDE computational errors
    if(isnan(h_M[0])){
        printf("\n\n PDE computation error. Possibly reduce the dT.\n\n");
        exit(1);
    }
    
    // Marking the start for the Profiler
    Profiler.Mark();
    
    // Clearing out the GridBin array
    for(int i=0;i<Grid_Size;i++) GridBin[i]=0;
    for(int i=0;i<MusselNo ;i++) MusselAlive[i]=ON;
    
    // Binning the mussel for each GridCell
    for(int i=0;i<MusselNo;i++)
    {
        GridX=min(Grid_Width -1,(int)(Mussel[i].x/ DomainSize * Grid_Width));
        GridY=min(Grid_Height-1,(int)(Mussel[i].y/ DomainSize * Grid_Height));
        
        Mussel[i].GridNr=GridY*Grid_Width+GridX;
        GridBin[Mussel[i].GridNr]++;
    }
    
    // Sorting the mussels per grid cell, needed when deleting mussels
    sort(Mussel.begin(), Mussel.end(), sortByGrid);
    
    // For each grid cell
    for(int j=0;j<Grid_Height;j++){
        for(int i=0;i<Grid_Width;i++){
            CurrentGridNr=j*Grid_Width+i;
            
            // How many mussels should be in the current cell, according to h_M
            MusselInCell=(int)(h_M[CurrentGridNr]*dX*dY*Mussels_Per_Gr+0.5);
            MusselChange = MusselInCell-GridBin[CurrentGridNr];
            if(MusselChange>0)
            {
                // Adding mussels to the current cell
                for (int k=0;k<MusselChange;k++)
                {
                    NewMussel.x=( (float)rand()/(float)RAND_MAX +(float)i )/Grid_Width  * DomainSize;
                    NewMussel.y=( (float)rand()/(float)RAND_MAX +(float)j )/Grid_Height * DomainSize;
                    NewMusselStack.push_back(NewMussel);
                }
            } else if(MusselChange<0)
            {
                // Removing "MusselChange" mussels from the field
                while(Mussel[counter].GridNr<CurrentGridNr) { counter++; }
                //Mussel.erase(Mussel.begin()+counter, Mussel.begin()+counter-MusselChange);
                for(int i=counter;i<counter-MusselChange;i++) MusselAlive[i]=OFF;
            }
        }
    }
    
    for(int i=0;i<MusselNo;i++)
    {
        if(MusselAlive[i]==ON)
        {
            OldMusselStack.push_back(Mussel[i]);
        }
    }
    
    OldMusselStack.insert(OldMusselStack.end(), NewMusselStack.begin(), NewMusselStack.end());
    Mussel=OldMusselStack;
    
    // Checking how many mussels should be in the field
    float M_Sum = 0;
    for(int i=0;i<Grid_Size;i++) M_Sum=M_Sum+h_M[i]*dX*dY;
    int TargetMusselNr = M_Sum*Mussels_Per_Gr;
    
    //Checking how many mussels are in the field
    int M_Nr = (int)Mussel.size();
    
    // If not enough mussels, mussels are added.
    if(M_Nr<TargetMusselNr){
        for (int k=0;k<(TargetMusselNr-M_Nr);k++){
            NewMussel.x=(float)rand()/(float)RAND_MAX*DomainSize;
            NewMussel.y=(float)rand()/(float)RAND_MAX*DomainSize;
            Mussel.push_back(NewMussel);
        }
    }
    
    // Updating the MusselNo variable
    MusselNo=(int)Mussel.size();
    
    // Registering the time spend
    Profiler.Register(2);
    
}

////////////////////////////////////////////////////////////////////////////////
// Translating the field of individials to gridded biomass
////////////////////////////////////////////////////////////////////////////////

void AdaptMusselGrid(float* &h_M, vector<Individual> Mussel, int MusselNo, ProfileClassCPU & Profiler)
{
    int GridX,GridY;
    int GridBin[Grid_Height*Grid_Width];
    
    // Marking the start for the Profiler
    Profiler.Mark();
    
    for(int i=0;i<Grid_Width*Grid_Height;i++) GridBin[i]=0;
    
    for(int i=0;i<MusselNo;i++)
    {
        GridX=min(Grid_Width-1,(int)(Mussel[i].x/DomainSize * Grid_Width));
        GridY=min(Grid_Height-1,(int)(Mussel[i].y/DomainSize * Grid_Height));
        
        Mussel[i].GridNr=GridY*Grid_Width+GridX;
        GridBin[Mussel[i].GridNr]++;
    }
    
    for(int i=0;i<Grid_Size;i++)
    {
        h_M[i] = (float)(GridBin[i]/Mussels_Per_Gr/dX/dY);
    }
    
    // Registering the time spend
    Profiler.Register(4);
    
}

