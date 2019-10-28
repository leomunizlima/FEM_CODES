#include "cylinder.h"

double CYLINDER_epresc(ParametersType *Parameters, double x, double y)
{
	double e, cv;
	
	cv = CYLINDER_cv(Parameters,x,y);

	e = cv + 0.5;
	
	e = 12.8664424218441;//7.64286; // rho * e = 12.8664424218441, para definir Mach = 0.38
	
	return e;
}

 
