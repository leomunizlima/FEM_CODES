# include "../01_CommonFiles/EulerEquations.h"

double NACA0512_rhopresc(double, double);
double NACA0512_v1presc(double, double);
double NACA0512_v2presc(double, double);
double NACA0512_epresc(ParametersType *, double, double);
double NACA0512_gamma(double, double);
double NACA0512_cv(ParametersType *, double, double);
double NACA0512_theta(double, double);
int NACA0512_InitialSolution(ParametersType *, NodeType *, double *);
void NACA0512_BC_no_penetrability(int, int, int, NodeType *, double, double, double, double, double, double, double, double, double, double, 
				  double [4][4], double [4][4], double [12][12], double [12][4], double [4][12], double [12], double [4], double [12], double [12], double [4], double [4]);



