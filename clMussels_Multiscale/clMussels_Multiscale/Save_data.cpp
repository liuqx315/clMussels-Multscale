//
//  HostFunctions.cpp
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
// Saving the PDE matrices, saves each stored array in a numbered file
////////////////////////////////////////////////////////////////////////////////

void Save_PDE_Matrices(char Prefix[], float* popA, float* popM, float* popS, int FileLab)
{
    // The location of the code is obtain from the __FILE__ macro
    const std::string SourcePath (__FILE__);
    const std::string PathName = SourcePath.substr (0,SourcePath.find_last_of("/")+1);
    
    float width_matrix = (float)Grid_Width;
    float height_matrix = (float)Grid_Height;
    
    char filename[2048];
    sprintf(filename,"%s%s%d.dat",PathName.c_str(),Prefix,FileLab);
    FILE* fp=fopen(filename,"w");
    
    fwrite(&width_matrix,sizeof(int),1,fp);
    fwrite(&height_matrix,sizeof(int),1,fp);
    
    fwrite(&popA[0],sizeof(float),Grid_Size,fp);
    fwrite(&popM[0],sizeof(float),Grid_Size,fp);
    fwrite(&popS[0],sizeof(float),Grid_Size,fp);
    
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
// Saving the X and Y vectors, saving at each stored set in a numbered file
////////////////////////////////////////////////////////////////////////////////

void Save_XY_Vectors(char Prefix[], std::vector<Individual> &Mussel, int len, int FileLab)
{
    // The location of the code is obtain from the __FILE__ macro
    const std::string SourcePath (__FILE__);
    const std::string PathName = SourcePath.substr (0,SourcePath.find_last_of("/")+1);
    
    char filename[2048];
    sprintf(filename,"%s%s%d.dat",PathName.c_str(),Prefix,FileLab);
    FILE* fp=fopen(filename,"w");
    
    float* Xc = (float*)malloc(sizeof(float)*len);
    float* Yc = (float*)malloc(sizeof(float)*len);
    
    for(int i=0;i<len;i++)
    {
        Xc[i]=Mussel[i].x;
        Yc[i]=Mussel[i].y;
    }
    
    int Domain_Size = DomainSize;
    int End_Time = EndTime;
    
    fwrite(&len,sizeof(int),1,fp);
    fwrite(&Domain_Size,sizeof(int),1,fp);
    fwrite(&End_Time, sizeof(int),1,fp);
    fwrite(&Xc[0],sizeof(float),len,fp);
    fwrite(&Yc[0],sizeof(float),len,fp);
    
    fclose(fp);
    free(Xc);
    free(Yc);
}
