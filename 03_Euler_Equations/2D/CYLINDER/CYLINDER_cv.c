#include "cylinder.h"

double CYLINDER_cv(ParametersType *Parameters, double x, double y)
{
	// calor específico
	double cv, M;
	M = Parameters->Mach;
	
	cv = 1.0/(M*M*0.56);
	
	//cv = 12.3664424218441;// 716.5 - Mach=2; //7.14286; // cv = 12.3664424218441 para definir mach = 0.38
	
	return cv;
	
}


