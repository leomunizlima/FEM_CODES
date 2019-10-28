#include "naca0012.h"

double NACA0012_epresc(ParametersType *Parameters, double x, double y)
{
	double e, cv;
	
	cv = NACA0012_cv(Parameters,x,y);

	e = cv + 0.5;
	
//	e = 1785714.7857143; 

	return e;
// rho * e = 2.28571->mach = 1.0 (Cv = 1.78571) %%% rho * e = 2.78795->mach = 0.8 (Cv = 2.28795) %%% rho * e = 7.642857->mach = 0.5 (Cv = 7.142857) %%% rho * e = 179.07143->mach = 0.1 (Cv=178.57143) %%% rho * e = 20.34127->mach = 0.3 (Cv=19.84127) %%% rho * e = 0.9464->mach = 2.0 (Cv = 716.5) %%% rho * e = 17857.64286->mach = 0.01 (Cv=17857.14286) %%% rho * e = 1785714.7857143->mach = 0.001 (Cv=1785714.2857143)		
}

 
