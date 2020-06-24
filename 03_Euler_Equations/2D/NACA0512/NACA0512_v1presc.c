#include "naca0512.h"

double NACA0512_v1presc(double x, double y){
	double v1;
	
	if ( fabs(x-0.0361133755718544)<=1e-15 && fabs(y+0.0307030292999)<=1e-15) {
		v1 = 0.0;
	}else{
		v1 = cos(0.261799387799149436538);
	}
	
	return v1;

}
 
