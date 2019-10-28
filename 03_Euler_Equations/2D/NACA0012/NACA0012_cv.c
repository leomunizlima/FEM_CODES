#include "naca0012.h"

double NACA0012_cv(ParametersType *Parameters, double x, double y)
{
	// calor específico
	double cv, M;
	M = Parameters->Mach;
	
	cv = 1.0/(M*M*0.56);
	
	//cv = 1785714.2857143; 
	
	return cv;
//Cv = 716.5->mach = 2.0 % Cv = 1.78571->mach = 1.0 % Cv = 2.28795->mach = 0.8 %  Cv = 7.142857->mach = 0.5 % Cv = 19.84127->mach = 0.3 % Cv = 178.57143->mach = 0.1 %  Cv = 17857.14286->mach = 0.01 % Cv = 1785714.2857143->mach = 0.001
	
	
	
}


