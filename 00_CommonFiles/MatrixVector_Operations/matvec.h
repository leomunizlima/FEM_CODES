#ifdef SSTranspEquation2D
	#include "../../01_SS_Transport_Equation/2D/01_CommonFiles/SSTranspEquation.h"
	int ebemv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int edemv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
#endif
#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"
	int ebemv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int edemv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
#endif
#ifdef EulerEquations2D
	#include "../../03_Euler_Equations/2D/01_CommonFiles/EulerEquations.h"
	int ebemvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int ebe2mvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int edemvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int ede2mvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
#endif
#ifdef FrozenEulerEquations2D
	#include "../../04_Frozen_Euler_Equations/2D/01_CommonFiles/FrozenEulerEquations.h"
#endif
#ifdef SSNavierStokesEquations2D
	#include "../../05_SS_NavierStokes_Equations/2D/01_CommonFiles/SSNavierStokesEquations.h"
#endif

int csrmv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);



