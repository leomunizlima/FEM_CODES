#include "EulerEquations.h"

int setStabilizationForm(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions,
						int (**Predictor)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *))
{
	
	if (strcasecmp(Parameters->StabilizationForm,"SUPG")==0){
		FemOtherFunctions->Build = Build_M_K_F_SUPG;
		if (strcasecmp(Parameters->TimeIntegration,"Predictor1")==0){		
			*Predictor = Predictor_Old;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor1_BDF")==0){		
			*Predictor = Predictor_Old_BDF;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor1_TRBDF2")==0){		
			*Predictor = Predictor_Old_TRBDF2;
		}
		else{
			printf("Time integration method is not defined correctly!\n");
			exit(1);
		}

		if (strcasecmp(Parameters->ShockCapture,"CAU")==0){
			FemFunctions->ShockCapture = Delta_CAU;
		}
		else if (strcasecmp(Parameters->ShockCapture,"YZBeta")==0){
			FemFunctions->ShockCapture = Delta_YZBeta;
		}	
		else{
			printf("Shock capture is not defined correctly!\n");
			exit(1);
		}
	}
	else if (strcasecmp(Parameters->StabilizationForm,"NMV1")==0||
		strcasecmp(Parameters->StabilizationForm,"NMV2")==0){
		if (strcasecmp(Parameters->TimeIntegration,"Predictor1")==0){		
			FemOtherFunctions->Build = Build_M_K_F_NMV_Transiente;
			*Predictor = Predictor_Old;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor1_BDF")==0){		
			FemOtherFunctions->Build = Build_M_K_F_NMV_Transiente;
			*Predictor = Predictor_Old_BDF;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor1_TRBDF2")==0){		
			FemOtherFunctions->Build = Build_M_K_F_NMV_Transiente;
			*Predictor = Predictor_Old_TRBDF2;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor2")==0){		
			FemOtherFunctions->Build = Build_M_F_NMV_Transiente;
			*Predictor = Predictor_New;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor2_BDF")==0){		
			FemOtherFunctions->Build = Build_M_F_NMV_Transiente;
			*Predictor = Predictor_New_BDF;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor2_TRBDF2")==0){		
			FemOtherFunctions->Build = Build_M_F_NMV_Transiente;
			*Predictor = Predictor_New_TRBDF2;
		}
		else{
			printf("Time integration method is not difined correctly!\n");
			exit(1);
		}
		
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
				FemFunctions->ShockCapture = Delta_YZBeta;
				FemFunctions->ShockCapture_Micro = Delta_YZBeta;
			}
			else if (strcasecmp(Parameters->StabilizationForm,"NMV2")==0)
			{
				FemFunctions->ShockCapture = Delta_YZBeta;
				FemFunctions->ShockCapture_Micro = Delta_YZBetaNMV;
			}

		}	
		else{
			printf("Shock capture is not defined correctly!\n");
			exit(1);
		}
	}
	else {
		printf("Stabilizantion form is not defined!\n");
		exit(1);
	}

	return 0;
}

