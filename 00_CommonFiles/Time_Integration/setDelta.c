#include "time_integration.h"

int setDelta(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	
	if (strcasecmp(Parameters->ShockCapture,"CAUNMV")==0){

		if (strcasecmp(Parameters->StabilizationForm,"NMV1")==0)
		{
			FemFunctions->ShockCapture = Delta_CAU_NMV;
			FemFunctions->ShockCapture_Micro = Delta_CAU_NMV;
		}
		else if (strcasecmp(Parameters->StabilizationForm,"NMV2")==0)
		{
			FemFunctions->ShockCapture = Delta_CAU_NMV;
			FemFunctions->ShockCapture_Micro = Delta_CAU_NMV;
		}
	}
	else if (strcasecmp(Parameters->ShockCapture,"YZBeta")==0){
	
		if (strcasecmp(Parameters->StabilizationForm,"NMV1")==0)
		{
			FemFunctions->ShockCapture = Delta_YZBetaNMV1;
			FemFunctions->ShockCapture_Micro = Delta_YZBetaNMV1;
		}
		else if (strcasecmp(Parameters->StabilizationForm,"NMV2")==0)
		{
			FemFunctions->ShockCapture = Delta_YZBetaNMV2;
			FemFunctions->ShockCapture_Micro = Delta_YZBetaNMV2_micro;
		}

	}	

	return 0;
}



