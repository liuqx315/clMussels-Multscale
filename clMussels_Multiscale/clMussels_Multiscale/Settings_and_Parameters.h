//
//  Settings_and_Parameters.h
//
//  Created by Johan Van de Koppel on 03-09-14.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

// Preprocessor directives

#ifndef SETTINGS_AND_PARAMETERS_H
#define SETTINGS_AND_PARAMETERS_H

#define ON                  1
#define OFF                 0

#define ProgressBarWidth    45

#define PDE_Device_No       1
#define IBM_Device_No       2

#define Profiling           OFF
#define CheckNumbers        OFF
#define CheckSorting        OFF
#define CheckHashing        OFF
#define CheckPRNG           OFF

// ---- Parameters determining the simulation setup ----------------------------
#define DomainSize          5000    // In centimeter
#define Mussels_Per_Gr      2.0     // 2-5 mussels per gram
#define CyclesPerTimestep   1250    // Number of mussel movement steps per 1 PDE timestep

#define Time      0             // 0      - Start time of the simulation
#define dT        0.0001        // 0.0001 - The timestep of the simulation
#define EndTime	  (24*180)      // 24*180 - hours     The time at which the simulation ends
#define NumFrames 200           // Number of times during the simulation that the data is stored
#define	MAX_STORE (NumFrames+1) //

#define ProgressBarWidth 45

#define Maximum_Nr_Elements 1<<22 // Maximal number of individuals the simulation can handle (power of 2)
#define Minimum_Nr_Elements 1<<14 // Minimal number of individuals the simulation should have (power of 2)

/* ---------- PDE GPU context definition -------------------------------------*/
#define WorkGroupSize       16    // Size of a workgroup in one dimension
#define GridScale           128   // Total number of gridcells in one dimension

// Thread block size
#define Block_Size_X        (WorkGroupSize)
#define Block_Size_Y        (WorkGroupSize)

// Number of blocks
/* I define the Block_Number_ensions of the matrix as product of two numbers
 Makes it easier to keep them a multiple of something (16, 32) when using CUDA*/
#define Block_Number_X      (GridScale/WorkGroupSize)
#define Block_Number_Y      (GridScale/WorkGroupSize)

// Matrix Block_Number_ensions
// (chosen as multiples of the thread block size for simplicity)
#define Grid_Width          (Block_Size_X * Block_Number_X)	// Matrix A width
#define Grid_Height         (Block_Size_Y * Block_Number_Y)	// Matrix A height
#define Grid_Size           (Grid_Width*Grid_Height)        	// Grid Size
#define Grid_Memory         (sizeof(float) * Grid_Size)     	// Memory needed

// ---- PDE model parameters ---------------------------------------------------

// Parameters		    Original value   Explanation and Units
#define EX	1000		//  1            - Speeding constant

#define D	0.0000		//  0.0005       - m2/h  The diffusion constant describing the movement of mussels
#define V	0.2*60*60	//  0.1*60*60    - The advection constant describing the flow of algae with the water

#define Aup	1.5			//  1.00         - g/m3     Algal concentration in upper layer; Oosterschelde data
#define hh	0.10		//  0.10         - m        Height of the lower layer; defined
#define ff	100			//  100          - m3/m3/h  eXchange rate with upper layer; Guestimated
#define	cc	0.1			//  0.1          - g/g/h    Maximal consumption rate of the mussels; Riisgard 2001
#define	ee	0.200		//  0.2          - g/g      Trophic efficiency of mussels; Lamie
#define	dM	0.004		//  0.02         - g/g/h    Density dependent mortality rate of the mussels; Calibrated
#define	kM	500;		//  150          - g/m2     Effect of density on mortality; Guestimated

#define PG	0.1
#define PKS	20.0

#define PK1	0.0001
#define PDS	0.0050
#define DS	0.0005

#define dX	(DomainSize/100.0/Grid_Width)    // 0.05
#define dY	(DomainSize/100.0/Grid_Height)   // 0.05

// ---- Parameters of the function relating density to movement speed ----------
#define P1  100     // Conversion, movement towards each other when too few mussels around
#define P2  -80     // Conversion, movement away from each other when too many mussels around
#define P3  3       // Speed when mussel is entirely alone

#define D1  2       // Size of the direct neighborhood
#define D2  6       // Size of the cluster neighborhood

// ---- Math functions etc ----------------------------------------------------
#define pi          3.14159265359
#define max(A,B)	(A>B?A:B)
#define min(A,B)	(A<B?A:B)

// ---- Struct definitions for the mussels and the Hash bins -------------------

// Name definitions
#define MUSSELS     101
#define ALGAE       102
#define SEDIMENT	102

#endif



