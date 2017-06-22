#include "../preconditioners.h"
/* p = A z
 * A = (D + w U)^{-1} (D + w L)^{-1}
 * (D + w L) * p = w
 * (D + w U) * w = z
 * z out
*/
int SSOR_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2; //auxiliar
	double z0, z1, z2; //auxiliar
	double w0, w1, w2; //auxiliar
	double p0, p1, p2; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **A = MatrixData->A; //matrix data
	double **inv = MatrixData->invone; //1.0 / (De + 1.0)

	double omega;
	omega = atof(&(Parameters->Preconditioner[3])); // In SORw, w is the weight coefficient

	for (I=0; I<neq; I++)
		z[I] = 0;

	z[neq] = 0;

	// OPTION 1: (D + wU) (D +wL)
	
	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];

		p0 = p[lm0];
		p1 = p[lm1];
		p2 = p[lm2];

		// (D + w L) * w = p
		w0 = p0 * inv[I][0];
		w1 = (p1 - omega * A[I][3] * w0) * inv[I][1];
		w2 = (p2 - omega * A[I][6] * w0 - omega * A[I][7] * w1) * inv[I][2];

		//(D + w U) * z = w
		z2 = w2 * inv[I][2];
		z1 = (w1 - omega * A[I][5] * z2) * inv[I][1];
		z0 = (w0 - omega * A[I][2] * z2 - omega * A[I][1] * z1) * inv[I][0];

		z[lm0] += z0;
		z[lm1] += z1;
		z[lm2] += z2;

		z[neq] = 0.0;
	}


	// OPTION 2: (D + wU) D^{-1} (D +wL)
	/*
	double x0, x1, x2;
	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
		lm1 = lm[I][1];
		lm2 = lm[I][2];

		p0 = p[lm0];
		p1 = p[lm1];
		p2 = p[lm2];

		// (D + w L) * w = p
		w0 = p0 * inv[I][0];
		w1 = (p1 - omega * A[I][3] * w0) * inv[I][1];
		w2 = (p2 - omega * A[I][6] * w0 - omega * A[I][7] * w1) * inv[I][2];

		//D^{-1} x = w
		x0 = A[I][0] * w0;
		x1 = A[I][4] * w1;
		x2 = A[I][8] * w2;

		//(D + w U) * z = w
		z2 = x2 * inv[I][2];
		z1 = (x1 - omega * A[I][5] * z2) * inv[I][1];
		z0 = (x0 - omega * A[I][2] * z2 - omega * A[I][1] * z1) * inv[I][0];

		z[lm0] += z0;
		z[lm1] += z1;
		z[lm2] += z2;

		z[neq] = 0.0;
	}
	*/

	return 0;
}
