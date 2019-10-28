#include "EulerEquations.h"

void setLocalPreconditioning(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Preconditioned,"VLR")==0)
		FemFunctions->Ax_Ay_calculations = VLR_Preconditioner_Ax_Ay_calculations;
	else if (strcasecmp(Parameters->Preconditioned,"WSCM")==0)
		FemFunctions->Ax_Ay_calculations = WSCM_Preconditioner_Ax_Ay_calculations;
	else if (strcasecmp(Parameters->Preconditioned,"NOT")==0)
		FemFunctions->Ax_Ay_calculations = NPreconditioned_Ax_Ay_calculations;
	else{
		printf("Local Preconditioning is not defined correctly!\n");
		exit(1);
	}	

}
