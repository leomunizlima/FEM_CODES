#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int SSOR_precond_EDE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int neq = Parameters->neq;
	int nedge = Parameters->nedge;
	
	if (tag==1){	
		MatrixData->invoneaux = mycalloc("invDeaux of 'SSOR_precond_EDE_setup'",2*nedge,sizeof(double));
		MatrixData->invone = mycalloc("invDe of 'SSOR_precond_EDE_setup'",nedge,sizeof(double*));
		for (I=0; I<nedge; I++)
			MatrixData->invone[I] = &(MatrixData->invoneaux[2*I]);
	}

	for (I=0; I<nedge; I++){
		MatrixData->invone[I][0] = 1.0/(1.0 + MatrixData->A[I][0]);
		MatrixData->invone[I][1] = 1.0/(1.0 + MatrixData->A[I][3]);
	}

	/* F preconditioning */
	double *faux = calloc((neq + 1), sizeof(double));
	for (I = 0; I < neq; I++){
		faux[I] = F[I];
	}
	faux[neq] = 0.0;
	F[neq] = 0.0;

	SSOR_precond_EDE (Parameters, MatrixData, FemStructs, faux, F);

	free(faux);

	return 0;

}

