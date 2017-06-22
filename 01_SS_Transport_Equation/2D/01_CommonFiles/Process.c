#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int i, neq = Parameters->neq;
	double NonLinearTolerance = Parameters->NonLinearTolerance; 
	double *u = FemStructs->u;
	double *uold, normdiff;	
	
	setProblem(Parameters,FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters, FemOtherFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions);
	
	uold = mycalloc("uold",sizeof(double),neq+1);	
	normdiff = NonLinearTolerance + 1;
	i = 0;
	while (normdiff > NonLinearTolerance){
		i++;	
		dcopy(neq, u, uold);
		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, i, FemStructs->F);
		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, FemStructs->F, FemStructs->u);

		daxpy(neq, -1, u, uold);
		normdiff = sqrt(ddot(neq,uold,uold));
		#ifdef debug
			printf("normdiff = %lf \n",normdiff);
		#endif
	}
	
	free(uold);

	return 0;
}


