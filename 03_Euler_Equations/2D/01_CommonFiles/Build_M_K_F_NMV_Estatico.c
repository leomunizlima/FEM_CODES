#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_M_K_F_NMV_Estatica(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int e, i, j, eNDOF;
	int J1, J2, J3, nel, neq;
	double x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21;
	double Area, twoArea, third = 1.0/3.0, sixth = 1.0/6.0, ninefortieth = 9.0 / 40.0;
	double Ub[4], Ue[12], dUb[4], dUe[12], *U, *dU;
	double delta, gamma, cv;
	double gradUx[4], gradUy[4];
	double Me[12][12], Ke[12][12], Ne[12][12], N[12], uBaux[4];
	double tolerance;
	double *R = FemStructs->F;
	double *uB_old = FemStructs->uB;
	double *delta_old = FemStructs->delta_old;
	double delta_t = Parameters->DeltaT_Build;
	double alpha = Parameters->Alpha_Build;
	int **lm = FemStructs->lm;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	tolerance = FemStructs->AuxBuild->tolerance;
	nel = Parameters->nel;
	neq = Parameters->neq;
	
	dzero(neq+1, R);
	setzeros(Parameters,MatrixData);

	U = (double*) mycalloc("U of 'Build_M_F_DD_Transiente'", 4*Parameters->nnodes, sizeof(double));
	dU = (double*) mycalloc("dU of 'Build_M_F_DD_Transiente'", 4*Parameters->nnodes, sizeof(double));
	eval_U_dU(Parameters,FemStructs,FemFunctions,U,dU);
		
	for (e = 0; e < nel; e++){
		// *** Global node that composes the element
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];

		// *** Nodal coordinates and differential operator
		x1 = Node[J1].x;
		x2 = Node[J2].x;
		x3 = Node[J3].x;
		
		y1 = Node[J1].y;
		y2 = Node[J2].y;
		y3 = Node[J3].y;
		
		gamma = FemFunctions->gamma(third*(x1+x2+x3), third*(y1+y2+y3));
		cv = FemFunctions->cv(Parameters, third*(x1+x2+x3), third*(y1+y2+y3));
		
		y23 = y2 - y3;
		y31 = y3 - y1;
		y12 = y1 - y2;

		x32 = x3 - x2;
		x13 = x1 - x3;
		x21 = x2 - x1;

		// *** Element Area
		twoArea = fabs(x21 * y31 - x13 * y12);
		Area = 0.5*twoArea; 


		// *** Fill Ue with the values ​​of the previous solution or workaround for each degree of liberty of the element node

		//Density
		Ue[0] = U[4*J1];
		Ue[4] = U[4*J2];
		Ue[8] = U[4*J3];

		dUe[0] = dU[4*J1];
		dUe[4] = dU[4*J2];
		dUe[8] = dU[4*J3];

		//x Velocity
		Ue[1] = U[4*J1+1];
		Ue[5] = U[4*J2+1];
		Ue[9] = U[4*J3+1];

		dUe[1] = dU[4*J1+1];
		dUe[5] = dU[4*J2+1];
		dUe[9] = dU[4*J3+1];

		//y Velocity
		Ue[2] = U[4*J1+2];
		Ue[6] = U[4*J2+2];
		Ue[10] = U[4*J3+2];

		dUe[2] = dU[4*J1+2];
		dUe[6] = dU[4*J2+2];
		dUe[10] = dU[4*J3+2];

		//Energy
		Ue[3] = U[4*J1+3];
		Ue[7] = U[4*J2+3];
		Ue[11] = U[4*J3+3];

		dUe[3] = dU[4*J1+3];
		dUe[7] = dU[4*J2+3];
		dUe[11] = dU[4*J3+3];


		// *** The triangle centroid
		for (i = 0; i < 4; i++){
			Ub[i] = (Ue[i] + Ue[i+4] + Ue[i+8]) * third;
			dUb[i] = (dUe[i] + dUe[i+4] + dUe[i+8]) * third;
		}

		// *** Fill uBaux with the values ​​of uB_old refers to the element
		eNDOF = e*NDOF;
		uBaux[0] = uB_old[eNDOF];
		uBaux[1] = uB_old[eNDOF + 1];
		uBaux[2] = uB_old[eNDOF + 2];
		uBaux[3] = uB_old[eNDOF + 3];
		
		// *** Calculation of Gradient (gradu = Bu)
		for (i = 0; i < 4; i++)
		{
			gradUx[i] = (Ue[i]*y23 + Ue[i+4]*y31 + Ue[i+8]*y12) / twoArea;
			gradUy[i] = (Ue[i]*x32 + Ue[i+4]*x13 + Ue[i+8]*x21) / twoArea;
		}
	
		//Maximum pressure variation on the element
		double p1, p2, p3, Dp12, Dp13, Dp23, Dpmax;
		
		p1 = (gamma - 1.0)*( Ue[3] - 1.0/(2*Ue[0])*(Ue[1]*Ue[1] + Ue[2]*Ue[2]) );		//presssure in node 1				
		p2 = (gamma - 1.0)*( Ue[7] - 1.0/(2*Ue[4])*(Ue[5]*Ue[5] + Ue[6]*Ue[6]) );		//presssure in node 2
		p3 = (gamma - 1.0)*( Ue[11] - 1.0/(2*Ue[8])*(Ue[9]*Ue[9] + Ue[10]*Ue[10]) );	//presssure in node 3	

		Dp12 = fabs(p1-p2);
		Dp13 = fabs(p1-p3);
		Dp23 = fabs(p2-p3);		
		
		if(Dp12 > Dp13)
			Dpmax = Dp12;
		else
			Dpmax = Dp13;
		if(Dpmax < Dp23)
			Dpmax = Dp23;

		// *** Jacobian matrices Ax and Ay calculations
		double Ax[4][4];
		double Ay[4][4];
	
		FemFunctions->Ax_Ay_calculations(cv,gamma,Parameters->Mach,Ub, Ax, Ay, Dpmax);

		// *** Inversa da Metrica Riemmaniana A0^(-1)
		double norma_U23 = (Ub[1] * Ub[1]) + (Ub[2] * Ub[2]);
		double rho_i = Ub[3] - ((norma_U23 * 0.5) / Ub[0]); 
		
		// Variaveis de entropia do baricentro
		double V2, V3, V4;
		V2 = Ub[1] / rho_i;
		V3 = Ub[2] / rho_i;
		V4 = (- Ub[0]) / rho_i;

		// A0^(-1) coefficients
		double rhoi_e    = - 1.0 / (rho_i * V4);
		double k1        = (V2*V2 + V3*V3) * 0.5 / V4;
		double A0[4][4];

		A0[0][0] = rhoi_e * (k1 * k1 + gamma);
		A0[0][1] = rhoi_e * (k1 * V2);
		A0[0][2] = rhoi_e * (k1 * V3);
		A0[0][3] = rhoi_e * ((k1 + 1.0) * V4);
		A0[1][0] = A0[0][1];
		A0[1][1] = rhoi_e * (V2 * V2 - V4);
		A0[1][2] = rhoi_e * (V2 * V3);
		A0[1][3] = rhoi_e * (V2 * V4);
		A0[2][0] = A0[0][2];
		A0[2][1] = A0[1][2];
		A0[2][2] = rhoi_e * (V3 * V3 - V4);
		A0[2][3] = rhoi_e * (V3 * V4);
		A0[3][0] = A0[0][3];
		A0[3][1] = A0[1][3];
		A0[3][2] = A0[2][3];
		A0[3][3] = rhoi_e * V4 * V4;      
		
		//  Delta do método DD
		delta = FemFunctions->ShockCapture(tolerance, delta_old, gradUx, gradUy, Ax, Ay, A0, dUb, y23, y31, y12, x32, x13, x21, twoArea, e, Parameters->invY, Ub); //COM Y^(-1) FIXO
			
		//------------------------------------------------------------------------------
		//  MONTAGENS DAS MATRIZES
		//------------------------------------------------------------------------------

		//Coeficientes utilizados na matriz Khh e KBh

		// Matriz Ka = y23Ax + x32Ay
		double Ka11, Ka12, Ka13, Ka14, Ka21, Ka22, Ka23, Ka24, Ka31, Ka32, Ka33, Ka34, Ka41, Ka42, Ka43, Ka44;
		Ka11 = y23 * Ax[0][0] + x32 * Ay[0][0];
		Ka12 = y23 * Ax[0][1] + x32 * Ay[0][1];
		Ka13 = y23 * Ax[0][2] + x32 * Ay[0][2];
		Ka14 = y23 * Ax[0][3] + x32 * Ay[0][3];
		Ka21 = y23 * Ax[1][0] + x32 * Ay[1][0];
		Ka22 = y23 * Ax[1][1] + x32 * Ay[1][1];
		Ka23 = y23 * Ax[1][2] + x32 * Ay[1][2];
		Ka24 = y23 * Ax[1][3] + x32 * Ay[1][3];
		Ka31 = y23 * Ax[2][0] + x32 * Ay[2][0];
		Ka32 = y23 * Ax[2][1] + x32 * Ay[2][1];
		Ka33 = y23 * Ax[2][2] + x32 * Ay[2][2];
		Ka34 = y23 * Ax[2][3] + x32 * Ay[2][3];
		Ka41 = y23 * Ax[3][0] + x32 * Ay[3][0];
		Ka42 = y23 * Ax[3][1] + x32 * Ay[3][1];
		Ka43 = y23 * Ax[3][2] + x32 * Ay[3][2];
		Ka44 = y23 * Ax[3][3] + x32 * Ay[3][3];

		// Matriz Kb = y31Ax + x13Ay
		double Kb11, Kb12, Kb13, Kb14, Kb21, Kb22, Kb23, Kb24, Kb31, Kb32, Kb33, Kb34, Kb41, Kb42, Kb43, Kb44;
		Kb11 = y31 * Ax[0][0] + x13 * Ay[0][0];
		Kb12 = y31 * Ax[0][1] + x13 * Ay[0][1];
		Kb13 = y31 * Ax[0][2] + x13 * Ay[0][2];
		Kb14 = y31 * Ax[0][3] + x13 * Ay[0][3];
		Kb21 = y31 * Ax[1][0] + x13 * Ay[1][0];
		Kb22 = y31 * Ax[1][1] + x13 * Ay[1][1];
		Kb23 = y31 * Ax[1][2] + x13 * Ay[1][2];
		Kb24 = y31 * Ax[1][3] + x13 * Ay[1][3];
		Kb31 = y31 * Ax[2][0] + x13 * Ay[2][0];
		Kb32 = y31 * Ax[2][1] + x13 * Ay[2][1];
		Kb33 = y31 * Ax[2][2] + x13 * Ay[2][2];
		Kb34 = y31 * Ax[2][3] + x13 * Ay[2][3];
		Kb41 = y31 * Ax[3][0] + x13 * Ay[3][0];
		Kb42 = y31 * Ax[3][1] + x13 * Ay[3][1];
		Kb43 = y31 * Ax[3][2] + x13 * Ay[3][2];
		Kb44 = y31 * Ax[3][3] + x13 * Ay[3][3];

		// Matriz Kc = y12Ax + x21Ay
		double Kc11, Kc12, Kc13, Kc14, Kc21, Kc22, Kc23, Kc24, Kc31, Kc32, Kc33, Kc34, Kc41, Kc42, Kc43, Kc44;
		Kc11 = y12 * Ax[0][0] + x21 * Ay[0][0];
		Kc12 = y12 * Ax[0][1] + x21 * Ay[0][1];
		Kc13 = y12 * Ax[0][2] + x21 * Ay[0][2];
		Kc14 = y12 * Ax[0][3] + x21 * Ay[0][3];
		Kc21 = y12 * Ax[1][0] + x21 * Ay[1][0];
		Kc22 = y12 * Ax[1][1] + x21 * Ay[1][1];
		Kc23 = y12 * Ax[1][2] + x21 * Ay[1][2];
		Kc24 = y12 * Ax[1][3] + x21 * Ay[1][3];
		Kc31 = y12 * Ax[2][0] + x21 * Ay[2][0];
		Kc32 = y12 * Ax[2][1] + x21 * Ay[2][1];
		Kc33 = y12 * Ax[2][2] + x21 * Ay[2][2];
		Kc34 = y12 * Ax[2][3] + x21 * Ay[2][3];
		Kc41 = y12 * Ax[3][0] + x21 * Ay[3][0];
		Kc42 = y12 * Ax[3][1] + x21 * Ay[3][1];
		Kc43 = y12 * Ax[3][2] + x21 * Ay[3][2];
		Kc44 = y12 * Ax[3][3] + x21 * Ay[3][3];

		// *** Matriz de Massa Mhh 12x12 que é igual a matriz de massa de Galerkin    
		double Mhh1, Mhh2;
		Mhh1 = Area * sixth;
		Mhh2 = Area / 12.0;
		
		// *** Matriz de Massa MhB 12x4 e MBh 4x12, onde os coeficientes são iguais, contudo atribuiremos zero a MhB (condensacao estatica ) e 
		// os devidos valores a MBh
		
		double MhB;
		MhB = 0.0;

		double MBh;	
		MBh = (3.0 * Area) / 20.0;
		
		// *** Matriz de Massa MBB 4x4
		double MBB;
		MBB = 0.0;
		
		// *** Matriz de Rigidez advinda de Galerkin Khhg 12x12 onde as 4 primeiras linhas repetem duas vezes, segue calculo das 4 primeiras linhas
		double Khhg11, Khhg12, Khhg13, Khhg14, Khhg21, Khhg22, Khhg23, Khhg24, Khhg31, Khhg32, Khhg33, Khhg34, Khhg41, Khhg42, Khhg43, Khhg44;
		Khhg11 = Ka11 * sixth;
		Khhg12 = Ka12 * sixth;
		Khhg13 = Ka13 * sixth;
		Khhg14 = Ka14 * sixth;
		Khhg21 = Ka21 * sixth;
		Khhg22 = Ka22 * sixth;
		Khhg23 = Ka23 * sixth;
		Khhg24 = Ka24 * sixth;
		Khhg31 = Ka31 * sixth;
		Khhg32 = Ka32 * sixth;
		Khhg33 = Ka33 * sixth;
		Khhg34 = Ka34 * sixth;
		Khhg41 = Ka41 * sixth;
		Khhg42 = Ka42 * sixth;
		Khhg43 = Ka43 * sixth;
		Khhg44 = Ka44 * sixth;
		
		double Khhg15, Khhg16, Khhg17, Khhg18, Khhg25, Khhg26, Khhg27, Khhg28, Khhg35, Khhg36, Khhg37, Khhg38, Khhg45, Khhg46, Khhg47, Khhg48;
		Khhg15 = Kb11 * sixth;
		Khhg16 = Kb12 * sixth;
		Khhg17 = Kb13 * sixth;
		Khhg18 = Kb14 * sixth;
		Khhg25 = Kb21 * sixth;
		Khhg26 = Kb22 * sixth;
		Khhg27 = Kb23 * sixth;
		Khhg28 = Kb24 * sixth;
		Khhg35 = Kb31 * sixth;
		Khhg36 = Kb32 * sixth;
		Khhg37 = Kb33 * sixth;
		Khhg38 = Kb34 * sixth;
		Khhg45 = Kb41 * sixth;
		Khhg46 = Kb42 * sixth;
		Khhg47 = Kb43 * sixth;
		Khhg48 = Kb44 * sixth;
		
		double Khhg19, Khhg110, Khhg111, Khhg112, Khhg29, Khhg210, Khhg211, Khhg212, Khhg39, Khhg310, Khhg311, Khhg312, Khhg49, Khhg410, Khhg411, Khhg412;
		Khhg19  = Kc11 * sixth;
		Khhg110 = Kc12 * sixth;
		Khhg111 = Kc13 * sixth;
		Khhg112 = Kc14 * sixth;
		Khhg29  = Kc21 * sixth;
		Khhg210 = Kc22 * sixth;
		Khhg211 = Kc23 * sixth;
		Khhg212 = Kc24 * sixth;
		Khhg39  = Kc31 * sixth;
		Khhg310 = Kc32 * sixth;
		Khhg311 = Kc33 * sixth;
		Khhg312 = Kc34 * sixth;
		Khhg49  = Kc41 * sixth;
		Khhg410 = Kc42 * sixth;
		Khhg411 = Kc43 * sixth;
		Khhg412 = Kc44 * sixth;		

		// *** Matriz de Rigidez advinda do DD Khhdd 12x12 onde a matriz é simétrica e bloco diagonal, segue calculo das 6 componentes
		double Khhdd11, Khhdd15, Khhdd19, Khhdd55, Khhdd59, Khhdd99;
		double delta4area = delta / (4.0 * Area);
		
		Khhdd11 = delta4area * (y23 * y23 + x32 * x32);
		Khhdd15 = delta4area * (y23 * y31 + x32 * x13);
		Khhdd19 = delta4area * (y23 * y12 + x32 * x21);
		Khhdd55 = delta4area * (y31 * y31 + x13 * x13);
		Khhdd59 = delta4area * (y31 * y12 + x13 * x21);
		Khhdd99 = delta4area * (y12 * y12 + x21 * x21);
		
		// *** Matriz de Rigidez KhB 12x4
		double a1 = ninefortieth * (   y23  + 2.0*y31 + 2.0*y12);
		double a2 = ninefortieth * (   x32  + 2.0*x13 + 2.0*x21);
		double b1 = ninefortieth * (2.0*y23 +   y31   + 2.0*y12);
		double b2 = ninefortieth * (2.0*x32 +   x13   + 2.0*x21);
		double c1 = ninefortieth * (2.0*y23 + 2.0*y31 +   y12  );
		double c2 = ninefortieth * (2.0*x32 + 2.0*x13 +   x21  );
		
		double KhB11, KhB12, KhB13, KhB14, KhB21, KhB22, KhB23, KhB24, KhB31, KhB32, KhB33, KhB34, KhB41, KhB42, KhB43, KhB44;
		
		KhB11 = a1*Ax[0][0] + a2*Ay[0][0];
		KhB12 = a1*Ax[0][1] + a2*Ay[0][1];
		KhB13 = a1*Ax[0][2] + a2*Ay[0][2];
		KhB14 = a1*Ax[0][3] + a2*Ay[0][3];
		KhB21 = a1*Ax[1][0] + a2*Ay[1][0];
		KhB22 = a1*Ax[1][1] + a2*Ay[1][1];		
		KhB23 = a1*Ax[1][2] + a2*Ay[1][2];
		KhB24 = a1*Ax[1][3] + a2*Ay[1][3];
		KhB31 = a1*Ax[2][0] + a2*Ay[2][0];
		KhB32 = a1*Ax[2][1] + a2*Ay[2][1];
		KhB33 = a1*Ax[2][2] + a2*Ay[2][2];
		KhB34 = a1*Ax[2][3] + a2*Ay[2][3];
		KhB41 = a1*Ax[3][0] + a2*Ay[3][0];
		KhB42 = a1*Ax[3][1] + a2*Ay[3][1];
		KhB43 = a1*Ax[3][2] + a2*Ay[3][2];
		KhB44 = a1*Ax[3][3] + a2*Ay[3][3];
		
		double KhB51, KhB52, KhB53, KhB54, KhB61, KhB62, KhB63, KhB64, KhB71, KhB72, KhB73, KhB74, KhB81, KhB82, KhB83, KhB84;
		
		KhB51 = b1*Ax[0][0] + b2*Ay[0][0];
		KhB52 = b1*Ax[0][1] + b2*Ay[0][1];
		KhB53 = b1*Ax[0][2] + b2*Ay[0][2];
		KhB54 = b1*Ax[0][3] + b2*Ay[0][3];
		KhB61 = b1*Ax[1][0] + b2*Ay[1][0];
		KhB62 = b1*Ax[1][1] + b2*Ay[1][1];		
		KhB63 = b1*Ax[1][2] + b2*Ay[1][2];
		KhB64 = b1*Ax[1][3] + b2*Ay[1][3];
		KhB71 = b1*Ax[2][0] + b2*Ay[2][0];
		KhB72 = b1*Ax[2][1] + b2*Ay[2][1];
		KhB73 = b1*Ax[2][2] + b2*Ay[2][2];
		KhB74 = b1*Ax[2][3] + b2*Ay[2][3];
		KhB81 = b1*Ax[3][0] + b2*Ay[3][0];
		KhB82 = b1*Ax[3][1] + b2*Ay[3][1];
		KhB83 = b1*Ax[3][2] + b2*Ay[3][2];
		KhB84 = b1*Ax[3][3] + b2*Ay[3][3];
		
		double KhB91, KhB92, KhB93, KhB94, KhB101, KhB102, KhB103, KhB104, KhB111, KhB112, KhB113, KhB114, KhB121, KhB122, KhB123, KhB124;
		
		KhB91 = c1*Ax[0][0] + c2*Ay[0][0];
		KhB92 = c1*Ax[0][1] + c2*Ay[0][1];
		KhB93 = c1*Ax[0][2] + c2*Ay[0][2];
		KhB94 = c1*Ax[0][3] + c2*Ay[0][3];
		KhB101 = c1*Ax[1][0] + c2*Ay[1][0];
		KhB102 = c1*Ax[1][1] + c2*Ay[1][1];		
		KhB103 = c1*Ax[1][2] + c2*Ay[1][2];
		KhB104 = c1*Ax[1][3] + c2*Ay[1][3];
		KhB111 = c1*Ax[2][0] + c2*Ay[2][0];
		KhB112 = c1*Ax[2][1] + c2*Ay[2][1];
		KhB113 = c1*Ax[2][2] + c2*Ay[2][2];
		KhB114 = c1*Ax[2][3] + c2*Ay[2][3];
		KhB121 = c1*Ax[3][0] + c2*Ay[3][0];
		KhB122 = c1*Ax[3][1] + c2*Ay[3][1];
		KhB123 = c1*Ax[3][2] + c2*Ay[3][2];
		KhB124 = c1*Ax[3][3] + c2*Ay[3][3];
		
		// *** Matriz de Rigidez KBh 4x12
		double KBh11, KBh12, KBh13, KBh14, KBh21, KBh22, KBh23, KBh24, KBh31, KBh32, KBh33, KBh34, KBh41, KBh42, KBh43, KBh44;
		KBh11 = ninefortieth * Ka11;
		KBh12 = ninefortieth * Ka12;
		KBh13 = ninefortieth * Ka13;
		KBh14 = ninefortieth * Ka14;
		KBh21 = ninefortieth * Ka21;
		KBh22 = ninefortieth * Ka22;
		KBh23 = ninefortieth * Ka23;
		KBh24 = ninefortieth * Ka24;
		KBh31 = ninefortieth * Ka31;
		KBh32 = ninefortieth * Ka32;
		KBh33 = ninefortieth * Ka33;
		KBh34 = ninefortieth * Ka34;
		KBh41 = ninefortieth * Ka41;
		KBh42 = ninefortieth * Ka42;
		KBh43 = ninefortieth * Ka43;
		KBh44 = ninefortieth * Ka44;
		
		double KBh15, KBh16, KBh17, KBh18, KBh25, KBh26, KBh27, KBh28, KBh35, KBh36, KBh37, KBh38, KBh45, KBh46, KBh47, KBh48;
		KBh15 = ninefortieth * Kb11;
		KBh16 = ninefortieth * Kb12;
		KBh17 = ninefortieth * Kb13;
		KBh18 = ninefortieth * Kb14;
		KBh25 = ninefortieth * Kb21;
		KBh26 = ninefortieth * Kb22;
		KBh27 = ninefortieth * Kb23;
		KBh28 = ninefortieth * Kb24;
		KBh35 = ninefortieth * Kb31;
		KBh36 = ninefortieth * Kb32;
		KBh37 = ninefortieth * Kb33;
		KBh38 = ninefortieth * Kb34;
		KBh45 = ninefortieth * Kb41;
		KBh46 = ninefortieth * Kb42;
		KBh47 = ninefortieth * Kb43;
		KBh48 = ninefortieth * Kb44;
		
		double KBh19, KBh110, KBh111, KBh112, KBh29, KBh210, KBh211, KBh212, KBh39, KBh310, KBh311, KBh312, KBh49, KBh410, KBh411, KBh412;
		KBh19 = ninefortieth * Kc11;
		KBh110 = ninefortieth * Kc12;
		KBh111 = ninefortieth * Kc13;
		KBh112 = ninefortieth * Kc14;
		KBh29 = ninefortieth * Kc21;
		KBh210 = ninefortieth * Kc22;
		KBh211 = ninefortieth * Kc23;
		KBh212 = ninefortieth * Kc24;
		KBh39 = ninefortieth * Kc31;
		KBh310 = ninefortieth * Kc32;
		KBh311 = ninefortieth * Kc33;
		KBh312 = ninefortieth * Kc34;
		KBh49 = ninefortieth * Kc41;
		KBh410 = ninefortieth * Kc42;
		KBh411 = ninefortieth * Kc43;
		KBh412 = ninefortieth * Kc44;
		
		// *** Matriz de Rigidez KBB 4x4
		double KBB;
		
		KBB = (81.0*delta*(y23*y23 + y31*y31 + y12*y12 + y23*y31 + y23*y12 + y31*y12 + x32*x32 + x13*x13 + x21*x21 + x32*x13 + x32*x21 + x13*x21))/(40.0*Area);
		
		// *** Coeficientes da Matriz de Massa [Me]
		// Me = Mhh - KhB*KBB^(-1)*MBh
		//
		// Calculando inversa de MBB + delta_t*KBB
		// double invBB = 1.0 / (MBB + delta_t*KBB); (mudando...)
		double invKBB = 1.0/KBB;
				
/*		double invBB__t__MhB_p_deltat_t_KhB11 = invBB * (MhB + delta_t*KhB11);
		double invBB__t__deltat_t_KhB12 = invBB * delta_t*KhB12;
		double invBB__t__deltat_t_KhB13 = invBB * delta_t*KhB13;
		double invBB__t__deltat_t_KhB14 = invBB * delta_t*KhB14;
		double invBB__t__deltat_t_KhB21 = invBB * delta_t*KhB21;
		double invBB__t__MhB_p_deltat_t_KhB22 = invBB * (MhB + delta_t*KhB22);
		double invBB__t__deltat_t_KhB23 = invBB * delta_t*KhB23;
		double invBB__t__deltat_t_KhB24 = invBB * delta_t*KhB24;
		double invBB__t__deltat_t_KhB31 = invBB * delta_t*KhB31;
		double invBB__t__deltat_t_KhB32 = invBB * delta_t*KhB32;
		double invBB__t__MhB_p_deltat_t_KhB33 = invBB * (MhB + delta_t*KhB33);
		double invBB__t__deltat_t_KhB34 = invBB * delta_t*KhB34;
		double invBB__t__deltat_t_KhB41 = invBB * delta_t*KhB41;
		double invBB__t__deltat_t_KhB42 = invBB * delta_t*KhB42;
		double invBB__t__deltat_t_KhB43 = invBB * delta_t*KhB43;
		double invBB__t__MhB_p_deltat_t_KhB44 = invBB * (MhB + delta_t*KhB44);
		double invBB__t__MhB_p_deltat_t_KhB51 = invBB * (MhB + delta_t*KhB51);
		double invBB__t__deltat_t_KhB52 = invBB * delta_t*KhB52;
		double invBB__t__deltat_t_KhB53 = invBB * delta_t*KhB53;
		double invBB__t__deltat_t_KhB54 = invBB * delta_t*KhB54;
		double invBB__t__deltat_t_KhB61 = invBB * delta_t*KhB61;
		double invBB__t__MhB_p_deltat_t_KhB62 = invBB * (MhB + delta_t*KhB62);
		double invBB__t__deltat_t_KhB63 = invBB * delta_t*KhB63;
		double invBB__t__deltat_t_KhB64 = invBB * delta_t*KhB64;
		double invBB__t__deltat_t_KhB71 = invBB * delta_t*KhB71;
		double invBB__t__deltat_t_KhB72 = invBB * delta_t*KhB72;
		double invBB__t__MhB_p_deltat_t_KhB73 = invBB * (MhB + delta_t*KhB73);
		double invBB__t__deltat_t_KhB74 = invBB * delta_t*KhB74;
		double invBB__t__deltat_t_KhB81 = invBB * delta_t*KhB81;
		double invBB__t__deltat_t_KhB82 = invBB * delta_t*KhB82;
		double invBB__t__deltat_t_KhB83 = invBB * delta_t*KhB83;
		double invBB__t__MhB_p_deltat_t_KhB84 = invBB * (MhB + delta_t*KhB84);
		double invBB__t__MhB_p_deltat_t_KhB91 = invBB * (MhB + delta_t*KhB91);
		double invBB__t__deltat_t_KhB92 = invBB * delta_t*KhB92;
		double invBB__t__deltat_t_KhB93 = invBB * delta_t*KhB93;
		double invBB__t__deltat_t_KhB94 = invBB * delta_t*KhB94;
		double invBB__t__deltat_t_KhB101 = invBB * delta_t*KhB101;
		double invBB__t__MhB_p_deltat_t_KhB102 = invBB * (MhB + delta_t*KhB102);
		double invBB__t__deltat_t_KhB103 = invBB * delta_t*KhB103;
		double invBB__t__deltat_t_KhB104 = invBB * delta_t*KhB104;
		double invBB__t__deltat_t_KhB111 = invBB * delta_t*KhB111;
		double invBB__t__deltat_t_KhB112 = invBB * delta_t*KhB112;
		double invBB__t__MhB_p_deltat_t_KhB113 = invBB * (MhB + delta_t*KhB113);
		double invBB__t__deltat_t_KhB114 = invBB * delta_t*KhB114;
		double invBB__t__deltat_t_KhB121 = invBB * delta_t*KhB121;
		double invBB__t__deltat_t_KhB122 = invBB * delta_t*KhB122;
		double invBB__t__deltat_t_KhB123 = invBB * delta_t*KhB123;
		double invBB__t__MhB_p_deltat_t_KhB124 = invBB * (MhB + delta_t*KhB124);
*/

		//no1 - no1
		Me[0][0] = Mhh1 - invKBB * KhB11 * MBh;
		Me[0][1] =      - invKBB * KhB12 * MBh;
		Me[0][2] =      - invKBB * KhB13 * MBh;
		Me[0][3] =      - invKBB * KhB14 * MBh;
		Me[1][0] =      - invKBB * KhB21 * MBh;
		Me[1][1] = Mhh1 - invKBB * KhB22 * MBh;
		Me[1][2] =      - invKBB * KhB23 * MBh;
		Me[1][3] =      - invKBB * KhB24 * MBh;
		Me[2][0] =      - invKBB * KhB31 * MBh;
		Me[2][1] =      - invKBB * KhB32 * MBh;
		Me[2][2] = Mhh1 - invKBB * KhB33 * MBh;
		Me[2][3] =      - invKBB * KhB34 * MBh;
		Me[3][0] =      - invKBB * KhB41 * MBh;
		Me[3][1] =      - invKBB * KhB42 * MBh;
		Me[3][2] =      - invKBB * KhB43 * MBh;
		Me[3][3] = Mhh1 - invKBB * KhB44 * MBh;

		//no1 - no2 
		Me[0][4] = Mhh2 - invKBB * KhB11 * MBh;
		Me[0][5] =      - invKBB * KhB12 * MBh;
		Me[0][6] =      - invKBB * KhB13 * MBh;
		Me[0][7] =      - invKBB * KhB14 * MBh;
		Me[1][4] =      - invKBB * KhB21 * MBh;
		Me[1][5] = Mhh2 - invKBB * KhB22 * MBh;
		Me[1][6] =      - invKBB * KhB23 * MBh;
		Me[1][7] =      - invKBB * KhB24 * MBh;
		Me[2][4] =      - invKBB * KhB31 * MBh;
		Me[2][5] =      - invKBB * KhB32 * MBh;
		Me[2][6] = Mhh2 - invKBB * KhB33 * MBh;
		Me[2][7] =      - invKBB * KhB34 * MBh;
		Me[3][4] =      - invKBB * KhB41 * MBh;
		Me[3][5] =      - invKBB * KhB42 * MBh;
		Me[3][6] =      - invKBB * KhB43 * MBh;
		Me[3][7] = Mhh2 - invKBB * KhB44 * MBh;

		//no1 - no3
		Me[0][8]  = Mhh2 - invKBB * KhB11 * MBh;
		Me[0][9]  =      - invKBB * KhB12 * MBh;
		Me[0][10] =      - invKBB * KhB13 * MBh;
		Me[0][11] =      - invKBB * KhB14 * MBh;
		Me[1][8]  =      - invKBB * KhB21 * MBh;
		Me[1][9]  = Mhh2 - invKBB * KhB22 * MBh;
		Me[1][10] =      - invKBB * KhB23 * MBh;
		Me[1][11] =      - invKBB * KhB24 * MBh;
		Me[2][8]  =      - invKBB * KhB31 * MBh;
		Me[2][9]  =      - invKBB * KhB32 * MBh;
		Me[2][10] = Mhh2 - invKBB * KhB33 * MBh;
		Me[2][11] =      - invKBB * KhB34 * MBh;
		Me[3][8]  =      - invKBB * KhB41 * MBh;
		Me[3][9]  =      - invKBB * KhB42 * MBh;
		Me[3][10] =      - invKBB * KhB43 * MBh;
		Me[3][11] = Mhh2 - invKBB * KhB44 * MBh;

		//no2 - no1
		Me[4][0] = Mhh2 - invKBB * KhB51 * MBh;
		Me[4][1] =      - invKBB * KhB52 * MBh;
		Me[4][2] =      - invKBB * KhB53 * MBh;
		Me[4][3] =      - invKBB * KhB54 * MBh;
		Me[5][0] =      - invKBB * KhB61 * MBh;
		Me[5][1] = Mhh2 - invKBB * KhB62 * MBh;
		Me[5][2] =      - invKBB * KhB63 * MBh;
		Me[5][3] =      - invKBB * KhB64 * MBh;
		Me[6][0] =      - invKBB * KhB71 * MBh;
		Me[6][1] =      - invKBB * KhB72 * MBh;
		Me[6][2] = Mhh2 - invKBB * KhB73 * MBh;
		Me[6][3] =      - invKBB * KhB74 * MBh;
		Me[7][0] =      - invKBB * KhB81 * MBh;
		Me[7][1] =      - invKBB * KhB82 * MBh;
		Me[7][2] =      - invKBB * KhB83 * MBh;
		Me[7][3] = Mhh2 - invKBB * KhB84 * MBh;

		//no2 - no2
		Me[4][4] = Mhh1 - invKBB * KhB51 * MBh;
		Me[4][5] =      - invKBB * KhB52 * MBh;
		Me[4][6] =      - invKBB * KhB53 * MBh;
		Me[4][7] =      - invKBB * KhB54 * MBh;
		Me[5][4] =      - invKBB * KhB61 * MBh;
		Me[5][5] = Mhh1 - invKBB * KhB62 * MBh;
		Me[5][6] =      - invKBB * KhB63 * MBh;
		Me[5][7] =      - invKBB * KhB64 * MBh;
		Me[6][4] =      - invKBB * KhB71 * MBh;
		Me[6][5] =      - invKBB * KhB72 * MBh;
		Me[6][6] = Mhh1 - invKBB * KhB73 * MBh;
		Me[6][7] =      - invKBB * KhB74 * MBh;
		Me[7][4] =      - invKBB * KhB81 * MBh;
		Me[7][5] =      - invKBB * KhB82 * MBh;
		Me[7][6] =      - invKBB * KhB83 * MBh;
		Me[7][7] = Mhh1 - invKBB * KhB84 * MBh;

		//no2 - no3
		Me[4][8]  = Mhh2 - invKBB * KhB51 * MBh;
		Me[4][9]  =      - invKBB * KhB52 * MBh;
		Me[4][10] =      - invKBB * KhB53 * MBh;
		Me[4][11] =      - invKBB * KhB54 * MBh;
		Me[5][8]  =      - invKBB * KhB61 * MBh;
		Me[5][9]  = Mhh2 - invKBB * KhB62 * MBh;
		Me[5][10] =      - invKBB * KhB63 * MBh;
		Me[5][11] =      - invKBB * KhB64 * MBh;
		Me[6][8]  =      - invKBB * KhB71 * MBh;
		Me[6][9]  =      - invKBB * KhB72 * MBh;
		Me[6][10] = Mhh2 - invKBB * KhB73 * MBh;
		Me[6][11] =      - invKBB * KhB74 * MBh;
		Me[7][8]  =      - invKBB * KhB81 * MBh;
		Me[7][9]  =      - invKBB * KhB82 * MBh;
		Me[7][10] =      - invKBB * KhB83 * MBh;
		Me[7][11] = Mhh2 - invKBB * KhB84 * MBh;

		//no3 - no1
		Me[8][0]  = Mhh2 - invKBB * KhB91  * MBh;
		Me[8][1]  =      - invKBB * KhB92  * MBh;
		Me[8][2]  =      - invKBB * KhB93  * MBh;
		Me[8][3]  =      - invKBB * KhB94  * MBh;
		Me[9][0]  =      - invKBB * KhB101 * MBh;
		Me[9][1]  = Mhh2 - invKBB * KhB102 * MBh;
		Me[9][2]  =      - invKBB * KhB103 * MBh;
		Me[9][3]  =      - invKBB * KhB104 * MBh;
		Me[10][0] =      - invKBB * KhB111 * MBh;
		Me[10][1] =      - invKBB * KhB112 * MBh;
		Me[10][2] = Mhh2 - invKBB * KhB113 * MBh;
		Me[10][3] =      - invKBB * KhB114 * MBh;
		Me[11][0] =      - invKBB * KhB121 * MBh;
		Me[11][1] =      - invKBB * KhB122 * MBh;
		Me[11][2] =      - invKBB * KhB123 * MBh;
		Me[11][3] = Mhh2 - invKBB * KhB124 * MBh;

		//no3 - no2
		Me[8][4]  = Mhh2 - invKBB * KhB91 * MBh;
		Me[8][5]  =      - invKBB * KhB92 * MBh;
		Me[8][6]  =      - invKBB * KhB93 * MBh;
		Me[8][7]  =      - invKBB * KhB94 * MBh;
		Me[9][4]  =      - invKBB * KhB101 * MBh;
		Me[9][5]  = Mhh2 - invKBB * KhB102 * MBh;
		Me[9][6]  =      - invKBB * KhB103 * MBh;
		Me[9][7]  =      - invKBB * KhB104 * MBh;
		Me[10][4] =      - invKBB * KhB111 * MBh;
		Me[10][5] =      - invKBB * KhB112 * MBh;
		Me[10][6] = Mhh2 - invKBB * KhB113 * MBh;
		Me[10][7] =      - invKBB * KhB114 * MBh;
		Me[11][4] =      - invKBB * KhB121 * MBh;
		Me[11][5] =      - invKBB * KhB122 * MBh;
		Me[11][6] =      - invKBB * KhB123 * MBh;
		Me[11][7] = Mhh2 - invKBB * KhB124 * MBh;

		//no3 - no3
		Me[8][8]   = Mhh1 - invKBB * KhB91 * MBh;
		Me[8][9]   =      - invKBB * KhB92 * MBh;
		Me[8][10]  =      - invKBB * KhB93 * MBh;
		Me[8][11]  =      - invKBB * KhB94 * MBh;
		Me[9][8]   =      - invKBB * KhB101 * MBh;
		Me[9][9]   = Mhh1 - invKBB * KhB102 * MBh;
		Me[9][10]  =      - invKBB * KhB103 * MBh;
		Me[9][11]  =      - invKBB * KhB104 * MBh;
		Me[10][8]  =      - invKBB * KhB111 * MBh;
		Me[10][9]  =      - invKBB * KhB112 * MBh;
		Me[10][10] = Mhh1 - invKBB * KhB113 * MBh;
		Me[10][11] =      - invKBB * KhB114 * MBh;
		Me[11][8]  =      - invKBB * KhB121 * MBh;
		Me[11][9]  =      - invKBB * KhB122 * MBh;
		Me[11][10] =      - invKBB * KhB123 * MBh;
		Me[11][11] = Mhh1 - invKBB * KhB124 * MBh;

		// *** Coeficientes da Matriz de Rigidez [Ke]
		// Ke = Khh - KhB*KBB^(-1)*KBh
		//no1 - no1
		Ke[0][0] = Khhg11 + Khhdd11 - invKBB * KhB11 * KBh11 - invKBB * KhB12 * KBh21 -  invKBB * KhB13 * KBh31    -    invKBB * KhB14 * KBh41;
		Ke[0][1] =      Khhg12      - invKBB * KhB11 * KBh12 -    invKBB * KhB12 * KBh22    -    invKBB * KhB13 * KBh32    -    invKBB * KhB14 * KBh42;
		Ke[0][2] =      Khhg13      - invKBB * KhB11 * KBh13 -    invKBB * KhB12 * KBh23    -    invKBB * KhB13 * KBh33    -    invKBB * KhB14 * KBh43;
		Ke[0][3] =      Khhg14      - invKBB * KhB11 * KBh14 -    invKBB * KhB12 * KBh24    -    invKBB * KhB13 * KBh34    -    invKBB * KhB14 * KBh44;
		Ke[1][0] =      Khhg21      -    invKBB * KhB21 * KBh11    - invKBB * KhB22 * KBh21 -    invKBB * KhB23 * KBh31    -    invKBB * KhB24 * KBh41;
		Ke[1][1] = Khhg22 + Khhdd11 -    invKBB * KhB21 * KBh12    - invKBB * KhB22 * KBh22 -    invKBB * KhB23 * KBh32    -    invKBB * KhB24 * KBh42;
		Ke[1][2] =      Khhg23      -    invKBB * KhB21 * KBh13    - invKBB * KhB22 * KBh23 -    invKBB * KhB23 * KBh33    -    invKBB * KhB24 * KBh43;
		Ke[1][3] =      Khhg24      -    invKBB * KhB21 * KBh14    - invKBB * KhB22 * KBh24 -    invKBB * KhB23 * KBh34    -    invKBB * KhB24 * KBh44;
		Ke[2][0] =      Khhg31      -    invKBB * KhB31 * KBh11    -    invKBB * KhB32 * KBh21    - invKBB * KhB33 * KBh31 -    invKBB * KhB34 * KBh41;
		Ke[2][1] =      Khhg32      -    invKBB * KhB31 * KBh12    -    invKBB * KhB32 * KBh22    - invKBB * KhB33 * KBh32 -    invKBB * KhB34 * KBh42; 
		Ke[2][2] = Khhg33 + Khhdd11 -    invKBB * KhB31 * KBh13    -    invKBB * KhB32 * KBh23    - invKBB * KhB33 * KBh33 -    invKBB * KhB34 * KBh43;
		Ke[2][3] =      Khhg34      -    invKBB * KhB31 * KBh14    -    invKBB * KhB32 * KBh24    - invKBB * KhB33 * KBh34 -    invKBB * KhB34 * KBh44;
		Ke[3][0] =      Khhg41      -    invKBB * KhB41 * KBh11    -    invKBB * KhB42 * KBh21    -    invKBB * KhB43 * KBh31    - invKBB * KhB44 * KBh41;
		Ke[3][1] =      Khhg42      -    invKBB * KhB41 * KBh12    -    invKBB * KhB42 * KBh22    -    invKBB * KhB43 * KBh32    - invKBB * KhB44 * KBh42; 
		Ke[3][2] =      Khhg43      -    invKBB * KhB41 * KBh13    -    invKBB * KhB42 * KBh23    -    invKBB * KhB43 * KBh33    - invKBB * KhB44 * KBh43; 
		Ke[3][3] = Khhg44 + Khhdd11 -    invKBB * KhB41 * KBh14    -    invKBB * KhB42 * KBh24    -    invKBB * KhB43 * KBh34    - invKBB * KhB44 * KBh44; 

		//no1 - no2
		Ke[0][4] = Khhg15 + Khhdd15 - invKBB * KhB11 * KBh15 -    invKBB * KhB12 * KBh25    -    invKBB * KhB13 * KBh35    -    invKBB * KhB14 * KBh45; 
		Ke[0][5] =      Khhg16      - invKBB * KhB11 * KBh16 -    invKBB * KhB12 * KBh26    -    invKBB * KhB13 * KBh36    -    invKBB * KhB14 * KBh46; 
		Ke[0][6] =      Khhg17      - invKBB * KhB11 * KBh17 -    invKBB * KhB12 * KBh27    -    invKBB * KhB13 * KBh37    -    invKBB * KhB14 * KBh47; 
		Ke[0][7] =      Khhg18      - invKBB * KhB11 * KBh18 -    invKBB * KhB12 * KBh28    -    invKBB * KhB13 * KBh38    -    invKBB * KhB14 * KBh48; 
		Ke[1][4] =      Khhg25      -    invKBB * KhB21 * KBh15    - invKBB * KhB22 * KBh25 -    invKBB * KhB23 * KBh35    -    invKBB * KhB24 * KBh45; 
		Ke[1][5] = Khhg26 + Khhdd15 -    invKBB * KhB21 * KBh16    - invKBB * KhB22 * KBh26 -    invKBB * KhB23 * KBh36    -    invKBB * KhB24 * KBh46; 
		Ke[1][6] =      Khhg27      -    invKBB * KhB21 * KBh17    - invKBB * KhB22 * KBh27 -    invKBB * KhB23 * KBh37    -    invKBB * KhB24 * KBh47; 
		Ke[1][7] =      Khhg28      -    invKBB * KhB21 * KBh18    - invKBB * KhB22 * KBh28 -    invKBB * KhB23 * KBh38    -    invKBB * KhB24 * KBh48; 
		Ke[2][4] =      Khhg35      -    invKBB * KhB31 * KBh15    -    invKBB * KhB32 * KBh25    - invKBB * KhB33 * KBh35 -    invKBB * KhB34 * KBh45; 
		Ke[2][5] =      Khhg36      -    invKBB * KhB31 * KBh16    -    invKBB * KhB32 * KBh26    - invKBB * KhB33 * KBh36 -    invKBB * KhB34 * KBh46; 
		Ke[2][6] = Khhg37 + Khhdd15 -    invKBB * KhB31 * KBh17    -    invKBB * KhB32 * KBh27    - invKBB * KhB33 * KBh37 -    invKBB * KhB34 * KBh47; 
		Ke[2][7] =      Khhg38      -    invKBB * KhB31 * KBh18    -    invKBB * KhB32 * KBh28    - invKBB * KhB33 * KBh38 -    invKBB * KhB34 * KBh48; 
		Ke[3][4] =      Khhg45      -    invKBB * KhB41 * KBh15    -    invKBB * KhB42 * KBh25    -    invKBB * KhB43 * KBh35    - invKBB * KhB44 * KBh45; 
		Ke[3][5] =      Khhg46      -    invKBB * KhB41 * KBh16    -    invKBB * KhB42 * KBh26    -    invKBB * KhB43 * KBh36    - invKBB * KhB44 * KBh46;  
		Ke[3][6] =      Khhg47      -    invKBB * KhB41 * KBh17    -    invKBB * KhB42 * KBh27    -    invKBB * KhB43 * KBh37    - invKBB * KhB44 * KBh47;  
		Ke[3][7] = Khhg48 + Khhdd15 -    invKBB * KhB41 * KBh18    -    invKBB * KhB42 * KBh28    -    invKBB * KhB43 * KBh38    - invKBB * KhB44 * KBh48; 

		//no1 - no3
		Ke[0][8]  = Khhg19 + Khhdd19  - invKBB * KhB11 * KBh19  -    invKBB * KhB12 * KBh29     -    invKBB * KhB13 * KBh39     -    invKBB * KhB14 * KBh49; 
		Ke[0][9]  =      Khhg110      - invKBB * KhB11 * KBh110 -    invKBB * KhB12 * KBh210    -    invKBB * KhB13 * KBh310    -    invKBB * KhB14 * KBh410;  
		Ke[0][10] =      Khhg111      - invKBB * KhB11 * KBh111 -    invKBB * KhB12 * KBh211    -    invKBB * KhB13 * KBh311    -    invKBB * KhB14 * KBh411; 
		Ke[0][11] =      Khhg112      - invKBB * KhB11 * KBh112 -    invKBB * KhB12 * KBh212    -    invKBB * KhB13 * KBh312    -    invKBB * KhB14 * KBh412; 
		Ke[1][8]  =      Khhg29       -    invKBB * KhB21 * KBh19     - invKBB * KhB22 * KBh29  -    invKBB * KhB23 * KBh39     -    invKBB * KhB24 * KBh49;  
		Ke[1][9]  = Khhg210 + Khhdd19 -    invKBB * KhB21 * KBh110    - invKBB * KhB22 * KBh210 -    invKBB * KhB23 * KBh310    -    invKBB * KhB24 * KBh410; 
		Ke[1][10] =      Khhg211      -    invKBB * KhB21 * KBh111    - invKBB * KhB22 * KBh211 -    invKBB * KhB23 * KBh311    -    invKBB * KhB24 * KBh411; 
		Ke[1][11] =      Khhg212      -    invKBB * KhB21 * KBh112    - invKBB * KhB22 * KBh212 -    invKBB * KhB23 * KBh312    -    invKBB * KhB24 * KBh412; 
		Ke[2][8]  =      Khhg39       -    invKBB * KhB31 * KBh19     -    invKBB * KhB32 * KBh29     - invKBB * KhB33 * KBh39  -    invKBB * KhB34 * KBh49;  
		Ke[2][9]  =      Khhg310      -    invKBB * KhB31 * KBh110    -    invKBB * KhB32 * KBh210    - invKBB * KhB33 * KBh310 -    invKBB * KhB34 * KBh410;  
		Ke[2][10] = Khhg311 + Khhdd19 -    invKBB * KhB31 * KBh111    -    invKBB * KhB32 * KBh211    - invKBB * KhB33 * KBh311 -    invKBB * KhB34 * KBh411;  
		Ke[2][11] =      Khhg312      -    invKBB * KhB31 * KBh112    -    invKBB * KhB32 * KBh212    - invKBB * KhB33 * KBh312 -    invKBB * KhB34 * KBh412; 
		Ke[3][8]  =      Khhg49       -    invKBB * KhB41 * KBh19     -    invKBB * KhB42 * KBh29     -    invKBB * KhB43 * KBh39     - invKBB * KhB44 * KBh49; 
		Ke[3][9]  =      Khhg410      -    invKBB * KhB41 * KBh110    -    invKBB * KhB42 * KBh210    -    invKBB * KhB43 * KBh310    - invKBB * KhB44 * KBh410; 
		Ke[3][10] =      Khhg411      -    invKBB * KhB41 * KBh111    -    invKBB * KhB42 * KBh211    -    invKBB * KhB43 * KBh311    - invKBB * KhB44 * KBh411;   
		Ke[3][11] = Khhg412 + Khhdd19 -    invKBB * KhB41 * KBh112    -    invKBB * KhB42 * KBh212    -    invKBB * KhB43 * KBh312    - invKBB * KhB44 * KBh412; 

		//no2 - no1
		Ke[4][0] = Khhg11 + Khhdd15 - invKBB * KhB51 * KBh11 -    invKBB * KhB52 * KBh21    -    invKBB * KhB53 * KBh31    -    invKBB * KhB54 * KBh41; 
		Ke[4][1] =      Khhg12      - invKBB * KhB51 * KBh12 -    invKBB * KhB52 * KBh22    -    invKBB * KhB53 * KBh32    -    invKBB * KhB54 * KBh42; 
		Ke[4][2] =      Khhg13      - invKBB * KhB51 * KBh13 -    invKBB * KhB52 * KBh23    -    invKBB * KhB53 * KBh33    -    invKBB * KhB54 * KBh43; 
		Ke[4][3] =      Khhg14      - invKBB * KhB51 * KBh14 -    invKBB * KhB52 * KBh24    -    invKBB * KhB53 * KBh34    -    invKBB * KhB54 * KBh44; 
		Ke[5][0] =      Khhg21      -    invKBB * KhB61 * KBh11    - invKBB * KhB62 * KBh21 -    invKBB * KhB63 * KBh31    -    invKBB * KhB64 * KBh41; 
		Ke[5][1] = Khhg22 + Khhdd15 -    invKBB * KhB61 * KBh12    - invKBB * KhB62 * KBh22 -    invKBB * KhB63 * KBh32    -    invKBB * KhB64 * KBh42; 
		Ke[5][2] =      Khhg23      -    invKBB * KhB61 * KBh13    - invKBB * KhB62 * KBh23 -    invKBB * KhB63 * KBh33    -    invKBB * KhB64 * KBh43; 
		Ke[5][3] =      Khhg24      -    invKBB * KhB61 * KBh14    - invKBB * KhB62 * KBh24 -    invKBB * KhB63 * KBh34    -    invKBB * KhB64 * KBh44; 
		Ke[6][0] =      Khhg31      -    invKBB * KhB71 * KBh11    -    invKBB * KhB72 * KBh21    - invKBB * KhB73 * KBh31 -    invKBB * KhB74 * KBh41; 
		Ke[6][1] =      Khhg32      -    invKBB * KhB71 * KBh12    -    invKBB * KhB72 * KBh22    - invKBB * KhB73 * KBh32 -    invKBB * KhB74 * KBh42; 
		Ke[6][2] = Khhg33 + Khhdd15 -    invKBB * KhB71 * KBh13    -    invKBB * KhB72 * KBh23    - invKBB * KhB73 * KBh33 -    invKBB * KhB74 * KBh43; 
		Ke[6][3] =      Khhg34      -    invKBB * KhB71 * KBh14    -    invKBB * KhB72 * KBh24    - invKBB * KhB73 * KBh34 -    invKBB * KhB74 * KBh44; 
		Ke[7][0] =      Khhg41      -    invKBB * KhB81 * KBh11    -    invKBB * KhB82 * KBh21    -    invKBB * KhB83 * KBh31    - invKBB * KhB84 * KBh41;
		Ke[7][1] =      Khhg42      -    invKBB * KhB81 * KBh12    -    invKBB * KhB82 * KBh22    -    invKBB * KhB83 * KBh32    - invKBB * KhB84 * KBh42; 
		Ke[7][2] =      Khhg43      -    invKBB * KhB81 * KBh13    -    invKBB * KhB82 * KBh23    -    invKBB * KhB83 * KBh33    - invKBB * KhB84 * KBh43; 
		Ke[7][3] = Khhg44 + Khhdd15 -    invKBB * KhB81 * KBh14    -    invKBB * KhB82 * KBh24    -    invKBB * KhB83 * KBh34    - invKBB * KhB84 * KBh44; 

		//no2 - no2
		Ke[4][4] = Khhg15 + Khhdd55 - invKBB * KhB51 * KBh15 -    invKBB * KhB52 * KBh25    -    invKBB * KhB53 * KBh35    -    invKBB * KhB54 * KBh45; 
		Ke[4][5] =      Khhg16      - invKBB * KhB51 * KBh16 -    invKBB * KhB52 * KBh26    -    invKBB * KhB53 * KBh36    -    invKBB * KhB54 * KBh46;  
		Ke[4][6] =      Khhg17      - invKBB * KhB51 * KBh17 -    invKBB * KhB52 * KBh27    -    invKBB * KhB53 * KBh37    -    invKBB * KhB54 * KBh47; 
		Ke[4][7] =      Khhg18      - invKBB * KhB51 * KBh18 -    invKBB * KhB52 * KBh28    -    invKBB * KhB53 * KBh38    -    invKBB * KhB54 * KBh48;  
		Ke[5][4] =      Khhg25      -    invKBB * KhB61 * KBh15    - invKBB * KhB62 * KBh25 -    invKBB * KhB63 * KBh35    -    invKBB * KhB64 * KBh45; 
		Ke[5][5] = Khhg26 + Khhdd55 -    invKBB * KhB61 * KBh16    - invKBB * KhB62 * KBh26 -    invKBB * KhB63 * KBh36    -    invKBB * KhB64 * KBh46;  
		Ke[5][6] =      Khhg27      -    invKBB * KhB61 * KBh17    - invKBB * KhB62 * KBh27 -    invKBB * KhB63 * KBh37    -    invKBB * KhB64 * KBh47; 
		Ke[5][7] =      Khhg28      -    invKBB * KhB61 * KBh18    - invKBB * KhB62 * KBh28 -    invKBB * KhB63 * KBh38    -    invKBB * KhB64 * KBh48; 
		Ke[6][4] =      Khhg35      -    invKBB * KhB71 * KBh15    -    invKBB * KhB72 * KBh25    - invKBB * KhB73 * KBh35 -    invKBB * KhB74 * KBh45; 
		Ke[6][5] =      Khhg36      -    invKBB * KhB71 * KBh16    -    invKBB * KhB72 * KBh26    - invKBB * KhB73 * KBh36 -    invKBB * KhB74 * KBh46;
		Ke[6][6] = Khhg37 + Khhdd55 -    invKBB * KhB71 * KBh17    -    invKBB * KhB72 * KBh27    - invKBB * KhB73 * KBh37 -    invKBB * KhB74 * KBh47; 
		Ke[6][7] =      Khhg38      -    invKBB * KhB71 * KBh18    -    invKBB * KhB72 * KBh28    - invKBB * KhB73 * KBh38 -    invKBB * KhB74 * KBh48;  
		Ke[7][4] =      Khhg45      -    invKBB * KhB81 * KBh15    -    invKBB * KhB82 * KBh25    -    invKBB * KhB83 * KBh35    - invKBB * KhB84 * KBh45; 
		Ke[7][5] =      Khhg46      -    invKBB * KhB81 * KBh16    -    invKBB * KhB82 * KBh26    -    invKBB * KhB83 * KBh36    - invKBB * KhB84 * KBh46; 
		Ke[7][6] =      Khhg47      -    invKBB * KhB81 * KBh17    -    invKBB * KhB82 * KBh27    -    invKBB * KhB83 * KBh37    - invKBB * KhB84 * KBh47;
		Ke[7][7] = Khhg48 + Khhdd55 -    invKBB * KhB81 * KBh18    -    invKBB * KhB82 * KBh28    -    invKBB * KhB83 * KBh38    - invKBB * KhB84 * KBh48;  

		//no2 - no3
		Ke[4][8]  = Khhg19 + Khhdd59  - invKBB * KhB51 * KBh19  -    invKBB * KhB52 * KBh29     -    invKBB * KhB53 * KBh39     -    invKBB * KhB54 * KBh49; 
		Ke[4][9]  =       Khhg110     - invKBB * KhB51 * KBh110 -    invKBB * KhB52 * KBh210    -    invKBB * KhB53 * KBh310    -    invKBB * KhB54 * KBh410; 
		Ke[4][10] =       Khhg111     - invKBB * KhB51 * KBh111 -    invKBB * KhB52 * KBh211    -    invKBB * KhB53 * KBh311    -    invKBB * KhB54 * KBh411; 
		Ke[4][11] =       Khhg112     - invKBB * KhB51 * KBh112 -    invKBB * KhB52 * KBh212    -    invKBB * KhB53 * KBh312    -    invKBB * KhB54 * KBh412; 
		Ke[5][8]  =       Khhg29      -    invKBB * KhB61 * KBh19     - invKBB * KhB62 * KBh29  -    invKBB * KhB63 * KBh39     -    invKBB * KhB64 * KBh49;  
		Ke[5][9]  = Khhg210 + Khhdd59 -    invKBB * KhB61 * KBh110    - invKBB * KhB62 * KBh210 -    invKBB * KhB63 * KBh310    -    invKBB * KhB64 * KBh410;  
		Ke[5][10] =       Khhg211     -    invKBB * KhB61 * KBh111    - invKBB * KhB62 * KBh211 -    invKBB * KhB63 * KBh311    -    invKBB * KhB64 * KBh411; 
		Ke[5][11] =       Khhg212     -    invKBB * KhB61 * KBh112    - invKBB * KhB62 * KBh212 -    invKBB * KhB63 * KBh312    -    invKBB * KhB64 * KBh412;
		Ke[6][8]  =       Khhg39      -    invKBB * KhB71 * KBh19     -    invKBB * KhB72 * KBh29     - invKBB * KhB73 * KBh39  -    invKBB * KhB74 * KBh49;  
		Ke[6][9]  =       Khhg310     -    invKBB * KhB71 * KBh110    -    invKBB * KhB72 * KBh210    - invKBB * KhB73 * KBh310 -    invKBB * KhB74 * KBh410; 
		Ke[6][10] = Khhg311 + Khhdd59 -    invKBB * KhB71 * KBh111    -    invKBB * KhB72 * KBh211    - invKBB * KhB73 * KBh311 -    invKBB * KhB74 * KBh411; 
		Ke[6][11] =       Khhg312     -    invKBB * KhB71 * KBh112    -    invKBB * KhB72 * KBh212    - invKBB * KhB73 * KBh312 -    invKBB * KhB74 * KBh412; 
		Ke[7][8]  =       Khhg49      -    invKBB * KhB81 * KBh19     -    invKBB * KhB82 * KBh29     -    invKBB * KhB83 * KBh39     - invKBB * KhB84 * KBh49; 
		Ke[7][9]  =       Khhg410     -    invKBB * KhB81 * KBh110    -    invKBB * KhB82 * KBh210    -    invKBB * KhB83 * KBh310    - invKBB * KhB84 * KBh410; 
		Ke[7][10] =       Khhg411     -    invKBB * KhB81 * KBh111    -    invKBB * KhB82 * KBh211    -    invKBB * KhB83 * KBh311    - invKBB * KhB84 * KBh411; 
		Ke[7][11] = Khhg412 + Khhdd59 -    invKBB * KhB81 * KBh112    -    invKBB * KhB82 * KBh212    -    invKBB * KhB83 * KBh312    - invKBB * KhB84 * KBh412; 
		
		//no3 - no1
		Ke[8][0]  = Khhg11 + Khhdd19 - invKBB * KhB91 * KBh11  -    invKBB * KhB92 * KBh21     -    invKBB * KhB93 * KBh31     -    invKBB * KhB94 * KBh41; 
		Ke[8][1]  =      Khhg12      - invKBB * KhB91 * KBh12  -    invKBB * KhB92 * KBh22     -    invKBB * KhB93 * KBh32     -    invKBB * KhB94 * KBh42; 
		Ke[8][2]  =      Khhg13      - invKBB * KhB91 * KBh13  -    invKBB * KhB92 * KBh23     -    invKBB * KhB93 * KBh33     -    invKBB * KhB94 * KBh43;
		Ke[8][3]  =      Khhg14      - invKBB * KhB91 * KBh14  -    invKBB * KhB92 * KBh24     -    invKBB * KhB93 * KBh34     -    invKBB * KhB94 * KBh44; 
		Ke[9][0]  =      Khhg21      -    invKBB * KhB101 * KBh11    - invKBB * KhB102 * KBh21 -    invKBB * KhB103 * KBh31    -    invKBB * KhB104 * KBh41; 
		Ke[9][1]  = Khhg22 + Khhdd19 -    invKBB * KhB101 * KBh12    - invKBB * KhB102 * KBh22 -    invKBB * KhB103 * KBh32    -    invKBB * KhB104 * KBh42;
		Ke[9][2]  =      Khhg23      -    invKBB * KhB101 * KBh13    - invKBB * KhB102 * KBh23 -    invKBB * KhB103 * KBh33    -    invKBB * KhB104 * KBh43; 
		Ke[9][3]  =      Khhg24      -    invKBB * KhB101 * KBh14    - invKBB * KhB102 * KBh24 -    invKBB * KhB103 * KBh34    -    invKBB * KhB104 * KBh44; 
		Ke[10][0] =      Khhg31      -    invKBB * KhB111 * KBh11    -    invKBB * KhB112 * KBh21    - invKBB * KhB113 * KBh31 -    invKBB * KhB114 * KBh41; 

		Ke[10][1] =      Khhg32      -    invKBB * KhB111 * KBh12    -    invKBB * KhB112 * KBh22    - invKBB * KhB113 * KBh32 -    invKBB * KhB114 * KBh42; 
		Ke[10][2] = Khhg33 + Khhdd19 -    invKBB * KhB111 * KBh13    -    invKBB * KhB112 * KBh23    - invKBB * KhB113 * KBh33 -    invKBB * KhB114 * KBh43; 
		Ke[10][3] =      Khhg34      -    invKBB * KhB111 * KBh14    -    invKBB * KhB112 * KBh24    - invKBB * KhB113 * KBh34 -    invKBB * KhB114 * KBh44; 
		Ke[11][0] =      Khhg41      -    invKBB * KhB121 * KBh11    -    invKBB * KhB122 * KBh21    -    invKBB * KhB123 * KBh31    - invKBB * KhB124 * KBh41; 
		Ke[11][1] =      Khhg42      -    invKBB * KhB121 * KBh12    -    invKBB * KhB122 * KBh22    -    invKBB * KhB123 * KBh32    - invKBB * KhB124 * KBh42; 
		Ke[11][2] =      Khhg43      -    invKBB * KhB121 * KBh13    -    invKBB * KhB122 * KBh23    -    invKBB * KhB123 * KBh33    - invKBB * KhB124 * KBh43; 
		Ke[11][3] = Khhg44 + Khhdd19 -    invKBB * KhB121 * KBh14    -    invKBB * KhB122 * KBh24    -    invKBB * KhB123 * KBh34    - invKBB * KhB124 * KBh44; 

		//no3 - no2
		Ke[8][4]  = Khhg15 + Khhdd59 - invKBB * KhB91 * KBh15  -    invKBB * KhB92 * KBh25     -    invKBB * KhB93 * KBh35     -    invKBB * KhB94 * KBh45; 
		Ke[8][5]  =      Khhg16      - invKBB * KhB91 * KBh16  -    invKBB * KhB92 * KBh26     -    invKBB * KhB93 * KBh36     -    invKBB * KhB94 * KBh46; 
		Ke[8][6]  =      Khhg17      - invKBB * KhB91 * KBh17  -    invKBB * KhB92 * KBh27     -    invKBB * KhB93 * KBh37     -    invKBB * KhB94 * KBh47; 
		Ke[8][7]  =      Khhg18      - invKBB * KhB91 * KBh18  -    invKBB * KhB92 * KBh28     -    invKBB * KhB93 * KBh38     -    invKBB * KhB94 * KBh48;  
		Ke[9][4]  =      Khhg25      -    invKBB * KhB101 * KBh15    - invKBB * KhB102 * KBh25 -    invKBB * KhB103 * KBh35    -    invKBB * KhB104 * KBh45; 
		Ke[9][5]  = Khhg26 + Khhdd59 -    invKBB * KhB101 * KBh16    - invKBB * KhB102 * KBh26 -    invKBB * KhB103 * KBh36    -    invKBB * KhB104 * KBh46; 
		Ke[9][6]  =      Khhg27      -    invKBB * KhB101 * KBh17    - invKBB * KhB102 * KBh27 -    invKBB * KhB103 * KBh37    -    invKBB * KhB104 * KBh47; 
		Ke[9][7]  =      Khhg28      -    invKBB * KhB101 * KBh18    - invKBB * KhB102 * KBh28 -    invKBB * KhB103 * KBh38    -    invKBB * KhB104 * KBh48; 
		Ke[10][4] =      Khhg35      -    invKBB * KhB111 * KBh15    -    invKBB * KhB112 * KBh25    - invKBB * KhB113 * KBh35 -    invKBB * KhB114 * KBh45; 
		Ke[10][5] =      Khhg36      -    invKBB * KhB111 * KBh16    -    invKBB * KhB112 * KBh26    - invKBB * KhB113 * KBh36 -    invKBB * KhB114 * KBh46;  
		Ke[10][6] = Khhg37 + Khhdd59 -    invKBB * KhB111 * KBh17    -    invKBB * KhB112 * KBh27    - invKBB * KhB113 * KBh37 -    invKBB * KhB114 * KBh47;  
		Ke[10][7] =      Khhg38      -    invKBB * KhB111 * KBh18    -    invKBB * KhB112 * KBh28    - invKBB * KhB113 * KBh38 -    invKBB * KhB114 * KBh48;  
		Ke[11][4] =      Khhg45      -    invKBB * KhB121 * KBh15    -    invKBB * KhB122 * KBh25    -    invKBB * KhB123 * KBh35    - invKBB * KhB124 * KBh45; 
		Ke[11][5] =      Khhg46      -    invKBB * KhB121 * KBh16    -    invKBB * KhB122 * KBh26    -    invKBB * KhB123 * KBh36    - invKBB * KhB124 * KBh46; 
		Ke[11][6] =      Khhg47      -    invKBB * KhB121 * KBh17    -    invKBB * KhB122 * KBh27    -    invKBB * KhB123 * KBh37    - invKBB * KhB124 * KBh47; 
		Ke[11][7] = Khhg48 + Khhdd59 -    invKBB * KhB121 * KBh18    -    invKBB * KhB122 * KBh28    -    invKBB * KhB123 * KBh38    - invKBB * KhB124 * KBh48; 

		//no3 - no3
		Ke[8][8]   = Khhg19 + Khhdd99  - invKBB * KhB91 * KBh19   -    invKBB * KhB92 * KBh29      -    invKBB * KhB93 * KBh39      -    invKBB * KhB94 * KBh49;  
		Ke[8][9]   =      Khhg110      - invKBB * KhB91 * KBh110  -    invKBB * KhB92 * KBh210     -    invKBB * KhB93 * KBh310     -    invKBB * KhB94 * KBh410; 
		Ke[8][10]  =      Khhg111      - invKBB * KhB91 * KBh111  -    invKBB * KhB92 * KBh211     -    invKBB * KhB93 * KBh311     -    invKBB * KhB94 * KBh411; 
		Ke[8][11]  =      Khhg112      - invKBB * KhB91 * KBh112  -    invKBB * KhB92 * KBh212     -    invKBB * KhB93 * KBh312     -    invKBB * KhB94 * KBh412;  
		Ke[9][8]   =      Khhg29       -    invKBB * KhB101 * KBh19     - invKBB * KhB102 * KBh29  -    invKBB * KhB103 * KBh39     -    invKBB * KhB104 * KBh49; 
		Ke[9][9]   = Khhg210 + Khhdd99 -    invKBB * KhB101 * KBh110    - invKBB * KhB102 * KBh210 -    invKBB * KhB103 * KBh310    -    invKBB * KhB104 * KBh410;  
		Ke[9][10]  =      Khhg211      -    invKBB * KhB101 * KBh111    - invKBB * KhB102 * KBh211 -    invKBB * KhB103 * KBh311    -    invKBB * KhB104 * KBh411; 
		Ke[9][11]  =      Khhg212      -    invKBB * KhB101 * KBh112    - invKBB * KhB102 * KBh212 -    invKBB * KhB103 * KBh312    -    invKBB * KhB104 * KBh412; 
		Ke[10][8]  =      Khhg39       -    invKBB * KhB111 * KBh19     -    invKBB * KhB112 * KBh29     - invKBB * KhB113 * KBh39  -    invKBB * KhB114 * KBh49; 
		Ke[10][9]  =      Khhg310      -    invKBB * KhB111 * KBh110    -    invKBB * KhB112 * KBh210    - invKBB * KhB113 * KBh310 -    invKBB * KhB114 * KBh410; 
		Ke[10][10] = Khhg311 + Khhdd99 -    invKBB * KhB111 * KBh111    -    invKBB * KhB112 * KBh211    - invKBB * KhB113 * KBh311 -    invKBB * KhB114 * KBh411; 
		Ke[10][11] =      Khhg312      -    invKBB * KhB111 * KBh112    -    invKBB * KhB112 * KBh212    - invKBB * KhB113 * KBh312 -    invKBB * KhB114 * KBh412;  
		Ke[11][8]  =      Khhg49       -    invKBB * KhB121 * KBh19     -    invKBB * KhB122 * KBh29     -    invKBB * KhB123 * KBh39     - invKBB * KhB124 * KBh49; 
		Ke[11][9]  =      Khhg410      -    invKBB * KhB121 * KBh110    -    invKBB * KhB122 * KBh210    -    invKBB * KhB123 * KBh310    - invKBB * KhB124 * KBh410; 
		Ke[11][10] =      Khhg411      -    invKBB * KhB121 * KBh111    -    invKBB * KhB122 * KBh211    -    invKBB * KhB123 * KBh311    - invKBB * KhB124 * KBh411; 
		Ke[11][11] = Khhg412 + Khhdd99 -    invKBB * KhB121 * KBh112    -    invKBB * KhB122 * KBh212    -    invKBB * KhB123 * KBh312    - invKBB * KhB124 * KBh412;

		// *** Coeficientes da Matriz [Ne]
		// Ne = (MhB - (MhB + delta_t*KhB)(MBB + delta_t*KBB)^(-1)MBB) / delta_t
		
		double invdelta_t = 1.0 / delta_t; 
		
		//no1
		Ne[0][0] = (MhB - invBB__t__MhB_p_deltat_t_KhB11 * MBB) * invdelta_t;
		Ne[0][1] = (    - invBB__t__deltat_t_KhB12 * MBB) * invdelta_t;
		Ne[0][2] = (    - invBB__t__deltat_t_KhB13 * MBB) * invdelta_t;
		Ne[0][3] = (    - invBB__t__deltat_t_KhB14 * MBB) * invdelta_t;
		Ne[1][0] = (    - invBB__t__deltat_t_KhB21 * MBB) * invdelta_t;
		Ne[1][1] = (MhB - invBB__t__MhB_p_deltat_t_KhB22 * MBB) * invdelta_t;
		Ne[1][2] = (    - invBB__t__deltat_t_KhB23 * MBB) * invdelta_t;
		Ne[1][3] = (    - invBB__t__deltat_t_KhB24 * MBB) * invdelta_t;
		Ne[2][0] = (    - invBB__t__deltat_t_KhB31 * MBB) * invdelta_t;
		Ne[2][1] = (    - invBB__t__deltat_t_KhB32 * MBB) * invdelta_t;
		Ne[2][2] = (MhB - invBB__t__MhB_p_deltat_t_KhB33 * MBB) * invdelta_t;
		Ne[2][3] = (    - invBB__t__deltat_t_KhB34 * MBB) * invdelta_t;
		Ne[3][0] = (    - invBB__t__deltat_t_KhB41 * MBB) * invdelta_t;
		Ne[3][1] = (    - invBB__t__deltat_t_KhB42 * MBB) * invdelta_t;
		Ne[3][2] = (    - invBB__t__deltat_t_KhB43 * MBB) * invdelta_t;
		Ne[3][3] = (MhB - invBB__t__MhB_p_deltat_t_KhB44 * MBB) * invdelta_t;

		//no2
		Ne[4][0] = (MhB - invBB__t__MhB_p_deltat_t_KhB51 * MBB) * invdelta_t;
		Ne[4][1] = (    - invBB__t__deltat_t_KhB52 * MBB) * invdelta_t;
		Ne[4][2] = (    - invBB__t__deltat_t_KhB53 * MBB) * invdelta_t;
		Ne[4][3] = (    - invBB__t__deltat_t_KhB54 * MBB) * invdelta_t;
		Ne[5][0] = (    - invBB__t__deltat_t_KhB61 * MBB) * invdelta_t;
		Ne[5][1] = (MhB - invBB__t__MhB_p_deltat_t_KhB62 * MBB) * invdelta_t;
		Ne[5][2] = (    - invBB__t__deltat_t_KhB63 * MBB) * invdelta_t;
		Ne[5][3] = (    - invBB__t__deltat_t_KhB64 * MBB) * invdelta_t;
		Ne[6][0] = (    - invBB__t__deltat_t_KhB71 * MBB) * invdelta_t;
		Ne[6][1] = (    - invBB__t__deltat_t_KhB72 * MBB) * invdelta_t;
		Ne[6][2] = (MhB - invBB__t__MhB_p_deltat_t_KhB73 * MBB) * invdelta_t;
		Ne[6][3] = (    - invBB__t__deltat_t_KhB74 * MBB) * invdelta_t;
		Ne[7][0] = (    - invBB__t__deltat_t_KhB81 * MBB) * invdelta_t;
		Ne[7][1] = (    - invBB__t__deltat_t_KhB82 * MBB) * invdelta_t;
		Ne[7][2] = (    - invBB__t__deltat_t_KhB83 * MBB) * invdelta_t;
		Ne[7][3] = (MhB - invBB__t__MhB_p_deltat_t_KhB84 * MBB) * invdelta_t;

		//no3
		Ne[8][0]  = (MhB - invBB__t__MhB_p_deltat_t_KhB91  * MBB) * invdelta_t;
		Ne[8][1]  = (    - invBB__t__deltat_t_KhB92  * MBB) * invdelta_t;
		Ne[8][2]  = (    - invBB__t__deltat_t_KhB93  * MBB) * invdelta_t;
		Ne[8][3]  = (    - invBB__t__deltat_t_KhB94  * MBB) * invdelta_t;
		Ne[9][0]  = (    - invBB__t__deltat_t_KhB101 * MBB) * invdelta_t;
		Ne[9][1]  = (MhB - invBB__t__MhB_p_deltat_t_KhB102 * MBB) * invdelta_t;
		Ne[9][2]  = (    - invBB__t__deltat_t_KhB103 * MBB) * invdelta_t;
		Ne[9][3]  = (    - invBB__t__deltat_t_KhB104 * MBB) * invdelta_t;
		Ne[10][0] = (    - invBB__t__deltat_t_KhB111 * MBB) * invdelta_t;
		Ne[10][1] = (    - invBB__t__deltat_t_KhB112 * MBB) * invdelta_t;
		Ne[10][2] = (MhB - invBB__t__MhB_p_deltat_t_KhB113 * MBB) * invdelta_t;
		Ne[10][3] = (    - invBB__t__deltat_t_KhB114 * MBB) * invdelta_t;
		Ne[11][0] = (    - invBB__t__deltat_t_KhB121 * MBB) * invdelta_t;
		Ne[11][1] = (    - invBB__t__deltat_t_KhB122 * MBB) * invdelta_t;
		Ne[11][2] = (    - invBB__t__deltat_t_KhB123 * MBB) * invdelta_t;
		Ne[11][3] = (MhB - invBB__t__MhB_p_deltat_t_KhB124 * MBB) * invdelta_t;


		//*****************************************
		double Re[12], MedUe[12], KeUe[12], Ae[12][12];

		for (i=0; i<12; i++){
			MedUe[i] = 0;
			KeUe[i] = 0;
			for (j=0; j<12; j++){
				MedUe[i] += Me[i][j]*dUe[j];
				KeUe[i] += Ke[i][j]*Ue[j];
				Ae[i][j] = Me[i][j] + alpha*delta_t*Ke[i][j];
			}
			Re[i] = - MedUe[i] - KeUe[i];
			N[i] = Ne[i][0]*uBaux[0] + Ne[i][1]*uBaux[1] + Ne[i][2]*uBaux[2] + Ne[i][3]*uBaux[3];
		}

		
		// Calculando novo uB
		uB_old[eNDOF] = invBB*(MBB*uBaux[0] - delta_t*MhB*(dUe[0] + dUe[4] + dUe[8]) - delta_t*(KBh11*Ue[0] + KBh12*Ue[1] + KBh13*Ue[2] + KBh14*Ue[3] + KBh15*Ue[4] +
										KBh16*Ue[5] + KBh17*Ue[6] + KBh18*Ue[7] + KBh19*Ue[8] + KBh110*Ue[9] + KBh111*Ue[10] + KBh112*Ue[11]));
		uB_old[eNDOF + 1] = invBB*(MBB*uBaux[1] - delta_t*MhB*(dUe[1] + dUe[5] + dUe[9]) - delta_t*(KBh21*Ue[0] + KBh22*Ue[1] + KBh23*Ue[2] + KBh24*Ue[3] + KBh25*Ue[4] + 
										KBh26*Ue[5] + KBh27*Ue[6] + KBh28*Ue[7] + KBh29*Ue[8] + KBh210*Ue[9] + KBh211*Ue[10] + KBh212*Ue[11]));
		uB_old[eNDOF + 2] = invBB*(MBB*uBaux[2] - delta_t*MhB*(dUe[2] + dUe[6] + dUe[10]) - delta_t*(KBh31*Ue[0] + KBh32*Ue[1] + KBh33*Ue[2] + KBh34*Ue[3] + KBh35*Ue[4] + 
										KBh36*Ue[5] + KBh37*Ue[6] + KBh38*Ue[7] + KBh39*Ue[8] + KBh310*Ue[9] + KBh311*Ue[10] + KBh312*Ue[11]));
		uB_old[eNDOF + 3] = invBB*(MBB*uBaux[3] - delta_t*MhB*(dUe[3] + dUe[7] + dUe[11]) - delta_t*(KBh41*Ue[0] + KBh42*Ue[1] + KBh43*Ue[2] + KBh44*Ue[3] + KBh45*Ue[4] + 
										KBh46*Ue[5] + KBh47*Ue[6] + KBh48*Ue[7] + KBh49*Ue[8] + KBh410*Ue[9] + KBh411*Ue[10] + KBh412*Ue[11]));

		// Assembly of global R
		for (i = 0; i < 12; i++)
			R[lm[e][i]] += N[i] + Re[i];
		R[neq] = 0;

		FemFunctions->assembly(Parameters, MatrixData, FemStructs, e, Ae);
		
	}//for elemento

	// Liberacao dos espacos alocados
	free(U);
	free(dU);

	
	return 0;
}// end build
