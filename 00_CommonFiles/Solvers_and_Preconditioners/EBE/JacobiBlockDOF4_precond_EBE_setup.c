#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int JacobiBlockDOF4_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int neq = Parameters->neq;
	int nel = Parameters->nel;
	double denominator;
	double **A = MatrixData->A;

	if (tag==1){
		MatrixData->invDeaux = mycalloc("invDeaux of 'Jacobi_precond_EBE_setup'",NDOF*NNOEL*nel,sizeof(double));
		MatrixData->invDe = mycalloc("invDe of 'Jacobi_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++)
			MatrixData->invDe[I] = &(MatrixData->invDeaux[NDOF*NNOEL*I]);
	}

	for (I = 0; I < nel; I++){
		denominator = (A[I][0] * (A[I][13] * (A[I][26] * A[I][39] - A[I][27] * A[I][38]) + A[I][14] * (A[I][27] * A[I][37] - A[I][25] * A[I][39]) + A[I][15] * (A[I][25] * A[I][38] - A[I][26] * A[I][37])) +
				A[I][1] * (A[I][12] * (A[I][27] * A[I][38] - A[I][26] * A[I][39]) + A[I][14] * (A[I][24] * A[I][39] - A[I][27] * A[I][36]) + A[I][15] * (A[I][26] * A[I][36] - A[I][24] * A[I][38])) +
				A[I][2] * (A[I][12] * (A[I][25] * A[I][39] - A[I][27] * A[I][37]) + A[I][13] * (A[I][27] * A[I][36] - A[I][24] * A[I][39]) + A[I][15] * (A[I][24] * A[I][37] - A[I][25] * A[I][36])) +
				A[I][3] * (A[I][12] * (A[I][26] * A[I][37] - A[I][25] * A[I][38]) + A[I][13] * (A[I][24] * A[I][38] - A[I][26] * A[I][36]) + A[I][14] * (A[I][25] * A[I][36] - A[I][24] * A[I][37])));

		MatrixData->invDe[I][0] = 1.0 / denominator;

		denominator = (A[I][52] * (A[I][65] * (A[I][78] * A[I][91] - A[I][79] * A[I][90]) + A[I][66] * (A[I][79] * A[I][89] - A[I][77] * A[I][91]) + A[I][67] * (A[I][77] * A[I][90] - A[I][78] * A[I][89])) +
				A[I][53] * (A[I][64] * (A[I][79] * A[I][90] - A[I][78] * A[I][91]) + A[I][66] * (A[I][76] * A[I][91] - A[I][79] * A[I][88]) + A[I][67] * (A[I][78] * A[I][88] - A[I][76] * A[I][90])) +
				A[I][54] * (A[I][64] * (A[I][77] * A[I][91] - A[I][79] * A[I][89]) + A[I][65] * (A[I][79] * A[I][88] - A[I][76] * A[I][91]) + A[I][67] * (A[I][76] * A[I][89] - A[I][77] * A[I][88])) +
				A[I][55] * (A[I][64] * (A[I][78] * A[I][89] - A[I][77] * A[I][90]) + A[I][65] * (A[I][76] * A[I][90] - A[I][78] * A[I][88]) + A[I][66] * (A[I][77] * A[I][88] - A[I][76] * A[I][89])));

		MatrixData->invDe[I][1] = 1.0 / denominator;

		denominator = (A[I][104] * (A[I][117] * (A[I][130] * A[I][143] - A[I][131] * A[I][142]) + A[I][118] * (A[I][131] * A[I][141] - A[I][129] * A[I][143]) + A[I][119] * (A[I][129] * A[I][142] - A[I][130] * A[I][141])) +
				A[I][105] * (A[I][116] * (A[I][131] * A[I][142] - A[I][130] * A[I][143]) + A[I][118] * (A[I][128] * A[I][143] - A[I][131] * A[I][140]) + A[I][119] * (A[I][130] * A[I][140] - A[I][128] * A[I][142])) +
				A[I][106] * (A[I][116] * (A[I][129] * A[I][143] - A[I][131] * A[I][141]) + A[I][117] * (A[I][131] * A[I][140] - A[I][128] * A[I][143]) + A[I][119] * (A[I][128] * A[I][141] - A[I][129] * A[I][140])) +
				A[I][107] * (A[I][116] * (A[I][130] * A[I][141] - A[I][129] * A[I][142]) + A[I][117] * (A[I][128] * A[I][142] - A[I][130] * A[I][140]) + A[I][118] * (A[I][129] * A[I][140] - A[I][128] * A[I][141])));

		MatrixData->invDe[I][2] = 1.0 / denominator;
	}
	/* F preconditioning */
	double *faux = calloc((neq + 1), sizeof(double));
	for (I = 0; I < neq; I++){
		faux[I] = F[I];
	}
	faux[neq] = 0.0;
	F[neq] = 0.0;

	JacobiBlockDOF4_precond_EBE (Parameters, MatrixData, FemStructs, faux, F);

	free(faux);

	return 0;
}
