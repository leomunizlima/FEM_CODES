# include "../01_CommonFiles/EulerEquations.h"

double SOD_rhopresc(double, double);
double SOD_v1presc(double, double);
double SOD_v2presc(double, double);
double SOD_epresc(ParametersType *, double, double);
double SOD_gamma(double, double);
double SOD_cv(ParametersType *, double, double);
int SOD_InitialSolution(ParametersType *, NodeType *, double *);
