#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "SSNavierStokesEquations.h"

int Postprocess(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	FILE *OutFile;
	char FileName[2000];
	

	/*************************************************************/
	//		Paraview output to file
	/************************************************************/	  
	Paraview_Output(Parameters, FemStructs, FemFunctions);

	/*************************************************************/


	/****************************************************************************************/
	// 			Printing final result		
	/****************************************************************************************/
	printf("\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	printf("Problem Title: %s_%s\n", Parameters->ProblemTitle, Parameters->Experiments);
	printf("Number of nodes: %d\n", Parameters->nnodes);
	printf("Number of elements: %d\n", Parameters->nel);
	printf("Number of equations: %d\n", Parameters->neq);
	printf("Stabilization form used: %s\n", Parameters->StabilizationForm);
	printf("Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		printf("Reordering: %s (bandwidth before: %d) (bandwidth after: %d)\n", Parameters->reordering, 
		Parameters->bandwidth_bef, Parameters->bandwidth_aft);
	}
	printf("Solver used: %s\n", Parameters->Solver);
	printf("Linear Solver tolerance used: %E\n", Parameters->SolverTolerance);
	printf("Nonlinear Solver tolerance used: %E\n", Parameters->SolverToleranceNonLin);
	printf("Preconditioner: %s\n", Parameters->Preconditioner);
	printf("Maximum number of solver iteration: %d\n", Parameters->LinearMaxIter);
	printf("Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
	printf("Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);
	printf("Number of nonlinear iterations: %d\n", Parameters-> NLiterations);
	printf("\n========================================================================\n\n");
	
	sprintf(FileName,"../03_output/%s_%s_%s_%s_%s_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->StabilizationForm,	Parameters->MatrixVectorProductScheme, Parameters->Preconditioner,Parameters->nnodes, Parameters->nel);
	
	OutFile = myfopen(FileName,"w");
	fprintf(OutFile, "\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	fprintf(OutFile, "Problem Title: %s_%s\n", Parameters->ProblemTitle, Parameters->Experiments);
	fprintf(OutFile, "Number of nodes: %d\n", Parameters->nnodes);
	fprintf(OutFile, "Number of elements: %d\n", Parameters->nel);
	fprintf(OutFile, "Number of equations: %d\n", Parameters->neq);
	fprintf(OutFile, "Stabilization form used: %s\n", Parameters->StabilizationForm);
	fprintf(OutFile, "Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		fprintf(OutFile,"Reordering: %s (bandwidth before: %d) (bandwidth after: %d)\n", Parameters->reordering, 
		Parameters->bandwidth_bef, Parameters->bandwidth_aft);
	}
	fprintf(OutFile, "Solver used: %s\n", Parameters->Solver);
	fprintf(OutFile, "Linear Solver tolerance used: %E\n", Parameters->SolverTolerance);
	fprintf(OutFile, "Nonlinear Solver tolerance used: %E\n", Parameters->SolverToleranceNonLin);
	fprintf(OutFile, "Maximum number of solver iteration: %d\n", Parameters->LinearMaxIter);
	fprintf(OutFile, "Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
	fprintf(OutFile, "Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);
	fprintf(OutFile, "Number of nonlinear iterations: %d\n", Parameters-> NLiterations);

	fprintf(OutFile, " \n Normas: |v_x| = %.4lf, |v_y| = %.4lf, |p| = %.4lf (in Postprocess.c)\n", Parameters->normu, Parameters->normv, Parameters->normp);
	fprintf(OutFile, "\n========================================================================\n\n");
	fclose(OutFile);
	/****************************************************************************************/
	
	
	/***************************************************************************************/
	//				Memory deallocation	
	/**************************************************************************************/
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE") == 0){		
		free(MatrixData->A);
		free(MatrixData->Aaux);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE") == 0){	
		int I, nel;

		nel = Parameters->nel;
		free(MatrixData->A);
		free(MatrixData->Aaux);
		for (I = 0; I < nel; I++){
			free(MatrixData->Scheme_by_Element[I]);
			free(MatrixData->order[I]);
		}
		free(MatrixData->Scheme_by_Element);
		free(MatrixData->order);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0){		
		int I, nel;

		nel = Parameters->nel;
		free(MatrixData->AA);
		free(MatrixData->JA);
		free(MatrixData->IA);
		for (I = 0; I < nel; I++){
			free(MatrixData->Scheme_by_Element[I]);
		}
		free(MatrixData->Diag);
		free(MatrixData->invDiag);
		free(MatrixData->Scheme_by_Element);
		free(MatrixData->order);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
		if ((strncmp(Parameters->Preconditioner,"ILU",3)==0)){
			SPARILU_clean(MatrixData->ILUp);
			SPARMAT_clean(MatrixData->mat);
			free(MatrixData->Ailu);
		}
	}
	
	
	free(MatrixData);
	free(FemStructs->Node);
	free(FemStructs->Element);
	free(FemStructs->F);
	free(FemStructs->u);
	free(FemStructs);
	free(FemFunctions);
	free(FemOtherFunctions);
	
	/***************************************************************************************/	
	
	return 0;
}



