#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int SSORDOF4_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
  int I;
  int neq = Parameters->neq;
  int nel = Parameters->nel;

  if (tag==1){
		MatrixData->invDeaux = mycalloc("invDeaux of 'Jacobi_precond_EBE_setup'",NDOF*NNOEL*nel,sizeof(double));
		MatrixData->invDe = mycalloc("invDe of 'Jacobi_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++)
			MatrixData->invDe[I] = &(MatrixData->invDeaux[NDOF*NNOEL*I]);
	}

  for (I = 0; I < nel; I++){
    MatrixData->invDe[I][0] = 1.0/(1.0+MatrixData->A[0][I][0]);
    MatrixData->invDe[I][1] = 1.0/(1.0+MatrixData->A[0][I][13]);
    MatrixData->invDe[I][2] = 1.0/(1.0+MatrixData->A[0][I][26]);
    MatrixData->invDe[I][3] = 1.0/(1.0+MatrixData->A[0][I][39]);

    MatrixData->invDe[I][4] = 1.0/(1.0+MatrixData->A[0][I][52]);
    MatrixData->invDe[I][5] = 1.0/(1.0+MatrixData->A[0][I][65]);
    MatrixData->invDe[I][6] = 1.0/(1.0+MatrixData->A[0][I][78]);
    MatrixData->invDe[I][7] = 1.0/(1.0+MatrixData->A[0][I][91]);

    MatrixData->invDe[I][8] = 1.0/(1.0+MatrixData->A[0][I][104]);
    MatrixData->invDe[I][9] = 1.0/(1.0+MatrixData->A[0][I][117]);
    MatrixData->invDe[I][10] = 1.0/(1.0+MatrixData->A[0][I][130]);
    MatrixData->invDe[I][11] = 1.0/(1.0+MatrixData->A[0][I][143]);
  }

  /* F preconditioning */
  double *faux = calloc((neq + 1), sizeof(double));
  for (I = 0; I < neq; I++){
    faux[I] = F[I];
  }
  faux[neq] = 0.0;
  F[neq] = 0.0;

  SSORDOF4_precond_EBE (Parameters, MatrixData, faux, F);

  free(faux);

  return 0;
}
