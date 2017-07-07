#include "scaling.h"
#include "preconditioners.h"

int NO_scaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	return 0;
}


int Left_scaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	int I, J, E, nel = Parameters->nel;
	int size = NNOEL*NDOF;
	int **lm = FemStructs->lm;
	double **A = MatrixData->A;
	double *invDiag = MatrixData->invDiag;	

	
	Diag_precond_EBE_setup(Parameters, MatrixData, FemStructs, 0, FemStructs->F);

	for (E=0; E<nel; E++){
		for (I=0; I<size; I++){
			for (J = 0; J<size; J++) 
				A[E][size*I+J] *= invDiag[lm[E][I]];	
		}
	}

	return 0;	
}


int LeftRight_scaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	int I, J, E, nel = Parameters->nel;
	int neq  = Parameters->neq;
	int size = NNOEL*NDOF;
	int **lm = FemStructs->lm;
	double **A = MatrixData->A;
	double *Diag = MatrixData->Diag;
	double *invDiag = MatrixData->invDiag;	
	double *F = FemStructs->F;
	
	Diag_precond_EBE_setup(Parameters, MatrixData, FemStructs, 0, F);

	for (E=0; E<nel; E++){
		for (I=0; I<size; I++){
			for (J = 0; J<size; J++) 
				A[E][size*I+J] *= sqrt(invDiag[lm[E][I]]*invDiag[lm[E][J]]);	
		}
	}

	for (I=0; I<neq; I++)
		F[I] *= sqrt(Diag[I]); //Da funcao anterior temos F divida pela diagonal Fi/dii , agora teremos Fi/sqrt(dii) pois (Fi/dii)*sqrt(di) = Fi/sqrt(dii). 
		

	return 0;
}

int NO_unscaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *x)
{
	return 0;
} 


int Left_unscaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *x)
{
	int I, neq = Parameters->neq; 
	double *invDiag = MatrixData->invDiag;

	
	for (I=0; I<neq; I++)
		x[I] *= sqrt(invDiag[I]);
		
	return 0;
}


