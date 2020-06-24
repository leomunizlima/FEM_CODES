#include "naca0512.h"

double NACA0512_v2presc(double x, double y){
	double v2;

	if ( fabs(x-0.0361133755718544)<=1e-15 && fabs(y+0.0307030292999)<=1e-15) {// || ( fabs(x - 1.0) <= 1e-4 && fabs(y) <= 1e-2 ) ){
		v2 = 0.0;
	}else{
		v2 = sin(0.261799387799149436538);
	}

	return v2;
} 


