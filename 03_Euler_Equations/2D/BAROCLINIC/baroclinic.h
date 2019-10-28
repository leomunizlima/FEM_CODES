# include "../01_CommonFiles/EulerEquations.h"

double BAROCLINIC_rhopresc(double, double);
double BAROCLINIC_v1presc(double, double);
double BAROCLINIC_v2presc(double, double);
double BAROCLINIC_epresc(ParametersType *,double, double);
double BAROCLINIC_gamma(double, double);
double BAROCLINIC_cv(ParametersType *, double, double);
int BAROCLINIC_InitialSolution(ParametersType *, NodeType *, double *);


