//
//  Profiler.cpp
//  clMussels_Multiscale
//
//  Created by Johan Van de Koppel on 03/01/2015.
//  Copyright (c) 2015 Johan Van de Koppel. All rights reserved.
//

#include <stdio.h>
#include <sys/time.h>
#include <vector>

#include "Profiler.h"

ProfileClassCPU::ProfileClassCPU(int Rec_Nr, int OnOff)
{
    
    for(int i=0;i<=Rec_Nr;i++)
    {
        Rec_Time.push_back(0.0);
    }
    
    IsOn = OnOff;
    
    gettimeofday(&Global_Time, NULL);
    gettimeofday(&Local_Time, NULL);
    
}

void ProfileClassCPU::Mark()
{
    //if(IsOn==1)
        gettimeofday(&Local_Time, NULL);
}

void ProfileClassCPU::Register(int Rec_Nr)
{
    gettimeofday(&Time_Now, NULL);
    double Time_Begin =Local_Time.tv_sec+(Local_Time.tv_usec/1000000.0);
    double Time_End   =Time_Now.tv_sec+(Time_Now.tv_usec/1000000.0);
    
    Rec_Time[Rec_Nr] += Time_End - Time_Begin;
}

double ProfileClassCPU::Report(int Rec_Nr)
{
    return Rec_Time[Rec_Nr];
}

double ProfileClassCPU::TimeSpanned()
{
    gettimeofday(&Time_Now, NULL);
    double Time_Begin =Global_Time.tv_sec+(Global_Time.tv_usec/1000000.0);
    double Time_End   =Time_Now.tv_sec+(Time_Now.tv_usec/1000000.0);
    
    double T = (Time_End - Time_Begin);
    
    return T;
}



