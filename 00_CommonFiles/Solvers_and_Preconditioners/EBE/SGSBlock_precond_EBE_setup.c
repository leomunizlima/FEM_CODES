#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int SGSBlock_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	SGSBlock_precond_EBE (Parameters, MatrixData, FemStructs, F, F);

	return 0;
}


