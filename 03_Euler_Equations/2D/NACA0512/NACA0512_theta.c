#include "naca0512.h"

double NACA0512_theta(double x, double y) 
{
	double theta;
	if( y > -1e-15 ) 
		theta = atan(0.0886749283154595/sqrt(x) - 0.0756000001025150 - 0.425687918985130*x + 0.520920909804897*x*x - 0.250184801216551*x*x*x);
	else
		theta = atan(-0.0886749283154595/sqrt(x) + 0.0756000001025150 + 0.425687918985130*x - 0.520920909804897*x*x + 0.250184801216551*x*x*x);
	/*
	if( y > -1e-10 ) 
		theta = atan(0.089070/sqrt(x) - 0.075600 - 0.421920*x + 0.511740*x*x - 0.243600*x*x*x);
	else
		theta = atan(-0.089070/sqrt(x) + 0.075600 + 0.421920*x - 0.511740*x*x + 0.243600*x*x*x);
	*/

	return theta;
}
/*
y = 0.594689181*(0.298222773*sqrt(x) - 0.127125232*x - 0.357907906*x*x +  0.291984971*x*x*x - 0.105174606*x*x*x*x);


y = 0.177349856630919*sqrt(x) - 0.0756000001025150*x - 0.212843959492565*x*x + 0.173640303268299*x*x*x - 0.0625462003041377*x*x*x*x


y' = 0.0886749283154595/sqrt(x) - 0.0756000001025150 - 0.425687918985130*x + 0.520920909804897*x*x - 0.250184801216551*x*x*x
*/
