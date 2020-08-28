#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

double Delta_NOT(ParametersType *Parameters,double *delta_old, double *gradUx, double *gradUy, double (*Ax)[4], double (*Ay)[4], double (*A0)[4], 
				double *dUb, double y23, double y31, double y12, double x32, double x13, double x21, double twoArea, int e, double *invY, double *Ub)
{
	
	return delta_old[e];
		
}

