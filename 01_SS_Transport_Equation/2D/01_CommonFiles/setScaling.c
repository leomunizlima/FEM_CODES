#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/scaling.h"

int setScaling(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Scaling,"NOT")==0){
		FemFunctions->scaling = NO_scaling;
		FemFunctions->unscaling = NO_unscaling;
	}
	else if (strcasecmp(Parameters->Scaling,"Left")==0){
		FemFunctions->scaling = Left_scaling;
		FemFunctions->unscaling = NO_unscaling;
	}
	else if (strcasecmp(Parameters->Scaling,"LeftRight")==0){
		FemFunctions->scaling = LeftRight_scaling;
		FemFunctions->unscaling = Left_unscaling;
	}
	else{
		printf("Scaling not defined correctly!\n");

	}


	return 0;
}
