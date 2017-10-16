#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int i, it, itmax;
	int iter, iterold, pricek;	
	double neq, delta1, delta2, epsilon, tol, norms, normu, etaold, eta0;	
	double *s, *Fold, *u, *F;

	
	setProblem(Parameters, FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters,FemOtherFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions);

	u = FemStructs->u;
	F = FemStructs->F;

	//=========Iteração de Newton===========
	neq = Parameters->neq;
	s = (double*) mycalloc("s of 'Process'", neq+1, sizeof(double));
	Fold = (double*) mycalloc("Fold of 'Process'", neq+1, sizeof(double));

	for(i=0; i<neq+1; i++){
		u[i] = 0.0;
		s[i] = 0.0;
		Fold[i] = 0.0;
	}

	//====== Constroi matriz e vetor força (Residuo) ======	
	FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
	delta1 = sqrt(ddot(neq, F, F));
	printf("\n Norma de F_0. |F_0| = %3.2E \n", delta1);
	it = 1;
	//====== Precondiona sistema ======
	FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, it, F);
	//precond_setup(Parameters, MatrixData, Node, Element, it, F, lm);

	delta1 = sqrt(ddot(neq, F, F));	
	delta2 = 1;	
	tol = 1e-3;	
	epsilon = tol*delta1;	
	printf("\n Norma de F_0. |F_0| = %3.2E === Saida de |F| = %3.2E  \n", delta1, epsilon);
	itmax = 1000;
	iter = 0;
	//Parameters->SolverTolerance = eta_0 (0.1, 0.5, 0.9) 
	eta0 = Parameters->SolverTolerance;
	
	while((delta1 > epsilon || delta2 > tol) && it < itmax){
		printf("===================================================");	
		printf("\n\n Eta: %3.2E, Iteracao Newton: %d, Pricek: %d \n\n", Parameters->SolverTolerance, it, pricek);
		it++;		
		iterold = iter;		
		dcopy(neq, F, Fold);   //Fould = F	
				
		//====== Resolve sitema linear ======		
		//Parameters->SolverTolerance = 1e-6;
		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, s);
		
		iter = Parameters->iterations;
		pricek = iter - iterold + 1; // tenho duvidas!!!!!
		daxpy(neq, 1, s, u);		//u = uold + s
				
		//====== Constroi matriz e vetor força ======		
		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
		//Build_K_F_SUPG_PSPG(Parameters, MatrixData, Node, Element, F, u, lm, assembly, end_assembly, ppresc, v1presc, v2presc);
		
		//====== Precondiona sistema ======
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, it, F);
		//precond_setup(Parameters, MatrixData, Node, Element, it, F, lm);

		//====== Atualiza eta ======
		etaold = Parameters->SolverTolerance;
		Parameters->SolverTolerance = eta_newton(Fold, F, etaold, pricek, it, epsilon, eta0, Parameters);
		
		//====== Condicoes de saida ======		
		delta1 = sqrt(ddot(neq, F, F));	
		norms = sqrt(ddot(neq, s, s));
		normu = sqrt(ddot(neq, u, u));
		delta2 = norms/normu;	
		printf("\n Norma de F. |F| = %3.2E  \n", delta1);		
		printf("\n         |s|/|u| = %3.2E  \n\n", delta2);
				
	}

	printf("===================================================\n");
	printf("\n Iterações de Newton = %d  \n", it-1);
		

	//SPARILU_clean(MatrixData->ILUp);	
	free(s);
	free(Fold);
	
	return 0;
}


















