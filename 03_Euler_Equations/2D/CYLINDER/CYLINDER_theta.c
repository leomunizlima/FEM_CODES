#include "cylinder.h"

double CYLINDER_theta(double x, double y)
{
	double theta;

	if (fabs(y)>=1e-12)
		theta = atan(-x/y); 
	else if (x<0)
		theta = 3.14159265359/2;
	else
		theta = -3.14159265359/2;


	//	printf("Math Error: division by zero!\n");
	//	exit(1);	

	return theta;
}



