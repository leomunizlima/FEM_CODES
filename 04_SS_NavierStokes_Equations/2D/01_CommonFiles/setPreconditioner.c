#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/preconditioners.h"

int setPreconditioner(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
        if (strcasecmp(Parameters->Preconditioner,"NOT")==0){
                FemFunctions->precond = NO_precond;
                FemFunctions->precondR = NO_precond;
                FemFunctions->precond_setup = NO_precond_setup;
/*                if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE2")==0){
                	*precond_setup = NOBlockDiag2_precond_EBE_setup;
		}
                if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE2")==0){
                	*precond_setup = NOBlockDiag2_precond_EDE_setup;
		}*/
        }
        else if (strcasecmp(Parameters->Preconditioner,"BlockDiag")==0){
                FemFunctions->precond = BlockDiagDOF3_precond;
                FemFunctions->precondR = NO_precond;
                if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
	                FemFunctions->precond_setup = BlockDiagDOF3_precond_EBE_setup;
                }
           /*     else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
	                *precond_setup = BlockDiag_precond_EDE_setup;
                }
                else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE2")==0){
	                *precond_setup = BlockDiag_precond_EBE_setup;
	                *precond = NO_precond;
		}		
                else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE2")==0){
	                *precond_setup = BlockDiag_precond_EDE_setup;
	                *precond = NO_precond;
		}		*/
		else{
			printf("Preconditioner definied only to EBE and EDE schemes\n");
			exit(1);
		}
			
	}
/*        else if (strcasecmp(Parameters->Preconditioner,"Diag")==0){
                *precond = Diag_precond;
                if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
	                *precond_setup = Diag_precond_EBE_setup;
                }
                else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
	                *precond_setup = Diag_precond_EDE_setup;
                }
		else{
			*precond_setup = Diag_precond_CSR_setup;
        	}
	}*/
        else if (strncmp(Parameters->Preconditioner,"ILU",3)==0){
                FemFunctions->precondR = NO_precond;
                if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
                	FemFunctions->precond = ILUp_precond;
			FemFunctions->precond_setup = ILUp_precond_setup;
		}
		else{
			printf("Preconditioner definied only to CSR scheme\n");
			exit(1);
		}
	}
	else {
                printf("Preconditioner not defined!\n");
                exit(1);

        }

        return 0;
}
