#include "../preconditioners.h"
/* Define preconditioner action on matrix-vector product */
/* p = A z
 * A = (LU)^{-1}
 * L * p = w
 * U * w = z
 * z out
*/
int LU_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z){
	int I;
	int lm0, lm1, lm2; //auxiliar
	double z0, z1, z2; //auxiliar
	double p0, p1, p2; //auxiliar
	double w0, w1, w2; //auxiliar
	int nel = Parameters->nel;  //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **LUe = MatrixData->LUe;

	//double **inv = MatrixData->invDe;
	//double **A = MatrixData->A;
	//double aux;

	for (I=0; I<neq; I++)
		z[I] = 0;

	z[neq] = 0;

	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];

		p0 = p[lm0];
		p1 = p[lm1];
		p2 = p[lm2];

		// L * p = w
		w0 = p0;
		w1 = (p1 - LUe[I][3] * w0);
		w2 = (p2 - LUe[I][6] * w0 - LUe[I][7] * w1);

		// U * w = z
		z2 = w2 * LUe[I][8];
		z1 = (w1 - LUe[I][5] * z2) / LUe[I][4];
		z0 = (w0 - LUe[I][2] * z2 - LUe[I][1] * z1) / LUe[I][0];

		z[lm0] += z0;
		z[lm1] += z1;
		z[lm2] += z2;

		z[neq] = 0.0;
	}

	return 0;
}
