#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int SGS_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int neq = Parameters->neq;
	int nel = Parameters->nel;

	if (tag==1){
		MatrixData->invoneaux = mycalloc("invDeaux of 'SGS_precond_EBE_setup'",NNOEL*nel,sizeof(double));
		MatrixData->invone = mycalloc("invDe of 'SGS_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++)
			MatrixData->invone[I] = &(MatrixData->invoneaux[NNOEL*I]);
	}

	for (I=0; I<nel; I++){
		MatrixData->invone[I][0] = 1.0/(1.0 + MatrixData->A[I][0]);
		MatrixData->invone[I][1] = 1.0/(1.0 + MatrixData->A[I][4]);
		MatrixData->invone[I][2] = 1.0/(1.0 + MatrixData->A[I][8]);
	}

	/* F preconditioning */
	double *faux = calloc((neq + 1), sizeof(double));
	for (I = 0; I < neq; I++){
		faux[I] = F[I];
	}
	faux[neq] = 0.0;
	F[neq] = 0.0;

	SGS_precond_EBE (Parameters, MatrixData, FemStructs, faux, F);

	free(faux);

	return 0;
}

