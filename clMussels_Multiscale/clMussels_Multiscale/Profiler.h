//
//  Profiler.h
//  clMussels_Multiscale
//
//  Created by Johan Van de Koppel on 03/01/2015.
//  Copyright (c) 2015 Johan Van de Koppel. All rights reserved.
//

#ifndef PROFILER_H
#define PROFILER_H

#include <stdio.h>
#include <sys/time.h>
#include <vector>

using namespace std;

class ProfileClassCPU{
    
private:
    
    vector<double> Rec_Time;
    struct timeval Global_Time;
    struct timeval Local_Time;
    struct timeval Time_Now;
    int IsOn;
    
public:
    ProfileClassCPU(int Rec_Nr, int OnOff);
    //~ProfileClassCPU();
    
    void Mark();
    void Register(int Rec_Nr);
    double Report(int Rec_Nr);
    double TimeSpanned();
    
};

#endif

