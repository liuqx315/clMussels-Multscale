#include "Settings_and_Parameters.h"

////////////////////////////////////////////////////////////////////////////////
// Laplacation operator definition, to calculate diffusive fluxes
////////////////////////////////////////////////////////////////////////////////

float LaplacianXY(__global float* pop, int row, int column)
{
	float retval;
	int current, left, right, top, bottom;
	float dx = dX;
	float dy = dY;
	
	current=row * Grid_Width + column;
	left=row * Grid_Width + column-1;
	right=row * Grid_Width + column+1;
	top=(row-1) * Grid_Width + column;
	bottom=(row+1) * Grid_Width + column;
    
	retval = ( (( pop[current] - pop[left] )/dx )
		      -(( pop[right]   - pop[current] )/dx )) / dx +
             ( (( pop[current] - pop[top] )/dy  )
              -(( pop[bottom]  - pop[current] )/dy ) ) / dy;
    
	return retval;
}

////////////////////////////////////////////////////////////////////////////////
// Gradient operator definition, to calculate advective fluxes
////////////////////////////////////////////////////////////////////////////////


float GradientY(__global float* pop, int row, int column)
{
	float retval;
	int current, top;
	float dy = dY;
	
	current=row * Grid_Width + column;
	top=(row-1) * Grid_Width + column;
	
	retval =  (( pop[current] - pop[top] )/dy );
    
	return retval;
}

////////////////////////////////////////////////////////////////////////////////
// Simulation kernel
////////////////////////////////////////////////////////////////////////////////

__kernel void MusselsKernel (__global float* A, __global float* M, __global float* S)
{
    
	float Dx=D*EX;
	float ex=ee*EX;
	float dMx=dM*EX;
    float PK1x=PK1*EX;
    float PDSx=PDS*EX;
    float DSx=DS*EX;
    float PkM = kM;
	
    size_t current = get_global_id(0);
    int row		= floor((float)current/(float)Grid_Width);
    int column	= current%Grid_Width;

	if (row > 0 && row < Grid_Width-1)
    {
        
        float dAdx	= GradientY(A, row, column);
        float d2Mdxy2 = LaplacianXY(M, row, column);
        float d2Sdxy2 = LaplacianXY(S, row, column);
        
        float Consumption = cc * A[current]*M[current]*(S[current]+PKS*PG)/(PKS+S[current]);
        
        A[current]=A[current]+(ff*(Aup-A[current]) - Consumption/hh - V*dAdx)*dT;
        M[current]=M[current]+(ex*Consumption - dMx*M[current] - Dx*d2Mdxy2)*dT;
        S[current]=S[current]+(PK1x*M[current]-PDSx*S[current]-DSx*d2Sdxy2)*dT;
        
    }
    
    //barrier(CLK_LOCAL_MEM_FENCE);
    
	// HANDLE Boundaries
    
	if(row==0)
		//do copy of first row = second last row
    {
        A[current]=A[(Grid_Height-2)*Grid_Width+column];
        M[current]=M[(Grid_Height-2)*Grid_Width+column];
        S[current]=S[(Grid_Height-2)*Grid_Width+column];
    }
    
	if(row==Grid_Height-1)
		//do copy of last row = second row
    {
        A[current]=A[1*Grid_Width+column];
        M[current]=M[1*Grid_Width+column];
        S[current]=S[1*Grid_Width+column];
    }
    
    if(column==0)
        //do copy of first row = second last row
    {
        A[current]=A[row * Grid_Width + Grid_Width-2];
        M[current]=M[row * Grid_Width + Grid_Width-2];
        S[current]=S[row * Grid_Width + Grid_Width-2];
    }
    
    if(column==Grid_Width-1)
        //do copy of last row = second row
    {
        A[current]=A[row * Grid_Width + 1];
        M[current]=M[row * Grid_Width + 1];
        S[current]=S[row * Grid_Width + 1];
    }
	
} // End Aridlandskernel

