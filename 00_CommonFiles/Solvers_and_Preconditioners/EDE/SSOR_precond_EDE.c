#include "../preconditioners.h"
/* z = A p
 * A = (D + U)^{-1} (D + L)^{-1}
 * (L + D) * p = w
 * (D + U) * z = p
 * z out
*/
int SSOR_precond_EDE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1; //auxiliar
	double z0, z1; //auxiliar
	double w0, w1; //auxiliar
	double p0, p1; //auxiliar
	int nedge = Parameters->nedge; //number of edges
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm2; //array LM
	double **A = MatrixData->A; //matrix data
//	double **inv = MatrixData->invone; //1.0 / (De + 1.0)

	double omega;
	omega = atof(&(Parameters->Preconditioner[4])); // In SORw, w is the weight coefficient

	for (I = 0; I < neq; I++){
		z[I] = 0.0;
	}

	for (I = 0; I < nedge; I++){
		lm0 = lm[I][0];
		lm1 = lm[I][1];

		p0 = p[lm0];
		p1 = p[lm1];

		/* (D + omega * L) * p = w */
		w0 = p0 / (1.0 + A[I][0]);
		w1 = (p1 - omega * A[I][2] * w0) / (1.0 + A[I][3]);

		/* (D + omega * U) * z = w */
		z1 = w1 / (1.0 + A[I][3]);
		z0 = (w0 - omega * A[I][1] * z1) / (1.0 + A[I][0]);

		z[lm0] += z0;
		z[lm1] += z1;

		z[neq] = 0.0;
	}

	return 0;
}
