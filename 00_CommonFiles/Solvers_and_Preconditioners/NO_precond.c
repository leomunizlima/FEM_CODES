#include "preconditioners.h"

int NO_precond(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int n = Parameters->neq;

	memcpy(z, p, n*sizeof(double));
		
	return 0;
}



