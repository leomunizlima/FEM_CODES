#ifdef SSTranspEquation2D
	#include "../../01_SS_Transport_Equation/2D/01_CommonFiles/SSTranspEquation.h"
	int NO_scaling(ParametersType *, MatrixDataType *, FemStructsType *); 
	int Left_scaling(ParametersType *, MatrixDataType *, FemStructsType *); 
	int LeftRight_scaling(ParametersType *, MatrixDataType *, FemStructsType *); 
	int NO_unscaling(ParametersType *, MatrixDataType *, FemStructsType *, double *); 
	int Left_unscaling(ParametersType *, MatrixDataType *, FemStructsType *, double *); 
#endif
#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"

#endif
