#include "../preconditioners.h"
/* z = A p
 * A = (D + U)^{-1} (D + L)^{-1}
 * (L + D) * p = w
 * (D + U) * z = p
 * z out
*/
int SGS_precond_EDE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1; //auxiliar
	double z0, z1; //auxiliar
	double w0, w1; //auxiliar
	double p0, p1; //auxiliar
	int nedge = Parameters->nedge; //number of edges
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm2; //array LM
	double **A = MatrixData->A;
	double **inv = MatrixData->invone; //1.0 / (De + 1.0)

	for (I = 0; I<neq; I++)
		z[I] = 0;
	z[neq] = 0;

	for (I = 0; I < nedge; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];

		p0 = p[lm0];
		p1 = p[lm1];

    		/* (D + L) * p = w */
		w0 = p0 * inv[I][0];
		w1 = (p1 - A[I][2] * w0) * inv[I][1];

		/* (D + U) * z = w */
		z1 = w1 / (1.0 + A[I][3]);
		z0 = (w0 - A[I][1] * z1) / (1.0 + A[I][0]);

		z[lm0] += z0;
		z[lm1] += z1;

		z[neq] = 0.0;
	}

	return 0;
}

