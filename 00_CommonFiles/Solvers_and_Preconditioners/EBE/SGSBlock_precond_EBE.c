#include "../preconditioners.h"

int SGSBlock_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11; //auxiliar
	double z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **A = MatrixData->A; //matrix data

	z[neq] = 0;

	for (I = nel-1; I >= 0; I--){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];
		lm3 = lm[I][3];
		lm4 = lm[I][4];
		lm5 = lm[I][5];
		lm6 = lm[I][6];
		lm7 = lm[I][7];
		lm8 = lm[I][8];
		lm9 = lm[I][9];
		lm10 = lm[I][10];
		lm11 = lm[I][11];

		z0 = z[lm0];
		z1 = z[lm1];
		z2 = z[lm2];
		z3 = z[lm3];
		z4 = z[lm4];
		z5 = z[lm5];
		z6 = z[lm6];
		z7 = z[lm7];
		z8 = z[lm8];
		z9 = z[lm9];
		z10 = z[lm10];
		z11 = z[lm11];

		z4 += - A[I][48]*z0;
		z5 += - A[I][61]*z1;
		z6 += - A[I][74]*z2;
		z7 += - A[I][74]*z3;

		z8 += - A[I][96]*z0 - A[I][100]*z4;
		z9 += - A[I][109]*z1 - A[I][113]*z5;
		z10 += - A[I][122]*z2 - A[I][126]*z6;
		z11 += - A[I][135]*z3 - A[I][139]*z7;

		z[lm4] = z4;
		z[lm5] = z5;
		z[lm6] = z6;
		z[lm7] = z7;
		z[lm8] = z8;
		z[lm9] = z9;
		z[lm10] = z10;
		z[lm11] = z11;

		z[neq] = 0.0;

	}


	return 0;
}

int SGSBlock_precondR_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11; //auxiliar
	double z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **A = MatrixData->A; //matrix data

	z[neq] = 0;


	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];
		lm3 = lm[I][3];
		lm4 = lm[I][4];
		lm5 = lm[I][5];
		lm6 = lm[I][6];
		lm7 = lm[I][7];
		lm8 = lm[I][8];
		lm9 = lm[I][9];
		lm10 = lm[I][10];
		lm11 = lm[I][11];

		z0 = z[lm0];
		z1 = z[lm1];
		z2 = z[lm2];
		z3 = z[lm3];
		z4 = z[lm4];
		z5 = z[lm5];
		z6 = z[lm6];
		z7 = z[lm7];
		z8 = z[lm8];
		z9 = z[lm9];
		z10 = z[lm10];
		z11 = z[lm11];

		z4 += -A[I][56]*z8;
		z5 += -A[I][69]*z9;
		z6 += -A[I][82]*z10;
		z7 += -A[I][95]*z11;
			
		z0 += -A[I][4]*z4 - A[I][8]*z8;
		z1 += -A[I][17]*z5 - A[I][21]*z9;
		z2 += -A[I][30]*z6 - A[I][34]*z10;
		z3 += -A[I][43]*z7 - A[I][47]*z11;
		z[lm0] = z0;
		z[lm1] = z1;
		z[lm2] = z2;
		z[lm3] = z3;
		z[lm4] = z4;
		z[lm5] = z5;
		z[lm6] = z6;
		z[lm7] = z7;

		z[neq] = 0.0;

	}

	return 0;
}


