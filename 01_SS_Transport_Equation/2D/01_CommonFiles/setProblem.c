#include "SSTranspEquation.h"
#include "../PUDIM/pudim.h"
#include "../CARTOLA/cartola.h"
#include "../TESTE/teste.h"
#include "../HEMKER/hemker.h"

int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{		
	if (strcasecmp(Parameters->ProblemTitle,"PUDIM")==0){
		FemFunctions->Condutivity = PUDIM_Condutivity;	
		FemFunctions->Font = PUDIM_Font;
		FemFunctions->Reaction = PUDIM_Reaction;
		FemFunctions->Velocity = PUDIM_Velocity;
		FemFunctions->upresc = PUDIM_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"CARTOLA")==0){
		FemFunctions->Condutivity = CARTOLA_Condutivity;	
		FemFunctions->Font = CARTOLA_Font;
		FemFunctions->Reaction = CARTOLA_Reaction;
		FemFunctions->Velocity = CARTOLA_Velocity;
		FemFunctions->upresc = CARTOLA_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"TESTE")==0){
		FemFunctions->Condutivity = TESTE_Condutivity;	
		FemFunctions->Font = TESTE_Font;
		FemFunctions->Reaction = TESTE_Reaction;
		FemFunctions->Velocity = TESTE_Velocity;
		FemFunctions->upresc = TESTE_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"HEMKER")==0){
		FemFunctions->Condutivity = HEMKER_Condutivity;	
		FemFunctions->Font = HEMKER_Font;
		FemFunctions->Reaction = HEMKER_Reaction;
		FemFunctions->Velocity = HEMKER_Velocity;
		FemFunctions->upresc = HEMKER_upresc;
	}
	else{
		printf("Problem not defined!\n");
		exit(1);
	}
	
	return 0;
}



