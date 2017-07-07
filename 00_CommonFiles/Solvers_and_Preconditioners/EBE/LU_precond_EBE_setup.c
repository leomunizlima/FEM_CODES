#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int LU_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int neq = Parameters->neq;
	int nel = Parameters->nel;
	double **LUe, **A;
	double m21, m31, m32, c22, c23, c32, c33, d33;

	if (tag == 1){
		MatrixData->LUeaux = mycalloc("Ueaux of 'LU_precond_EBE_setup'", 9*nel,sizeof(double));
		MatrixData->LUe = mycalloc("Ue of 'LU_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++){
					MatrixData->LUe[I] = &(MatrixData->LUeaux[9*I]);
		}
	}
	A = MatrixData->A;
	LUe = MatrixData->LUe;

	for (I = 0; I < nel; I++){
		m21 = -A[I][3]/A[I][0];
		m31 = -A[I][6]/A[I][0];
		c22 = A[I][4] + m21*A[I][1];
		c23	= A[I][5] + m21*A[I][2];
		c32 = A[I][7] + m31*A[I][1];
		c33 = A[I][8] + m31*A[I][2];
		m32 = -c32/c22;
		d33 = c33 + m32*c23;
		LUe[I][0] = A[I][0];
		LUe[I][1] = A[I][1];
		LUe[I][2] = A[I][2];
		LUe[I][3] = m21;
		LUe[I][4] = c22;
		LUe[I][5] = c23;
		LUe[I][6] = m31;
		LUe[I][7] = m32;
		LUe[I][8] = d33;
	}

	/* F preconditioning */
	double *faux = calloc((neq + 1), sizeof(double));
	for (I = 0; I < neq; I++){
		faux[I] = F[I];
	}
	faux[neq] = 0.0;
	F[neq] = 0.0;

	LU_precond_EBE (Parameters, MatrixData, FemStructs, faux, F);

	free(faux);

	return 0;

}
