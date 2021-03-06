#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"


int Build_M_F_NMV_Transiente(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int e, i, j, k, neq, nel;
	int J1, J2, J3;
	double x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21;
	double Area, twoArea, third = 1.0/3.0, sixth = 1.0/6.0, ninefortieth = 9.0 / 40.0;
	double Ub[4], dUb[4];
	double delta, deltaNMV, gamma, cv;
	double gradUx[4], gradUy[4];
	double **M2_out, **R1_out, **R2_out, *invN2_out, *U, *dU, *uB, *duB;
	double Me[12][12], Ue[12], dUe[12];
	double *delta_old = FemStructs->delta_old;
	double *deltaNMV_old = FemStructs->deltaNMV_old;
	double *F = FemStructs->F;
	double delta_t = Parameters->DeltaT_Build;
	double alpha = Parameters->Alpha_Build;
	int **lm = FemStructs->lm;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	//double delta_soma=0.0, delta_max=0.001;

	neq = Parameters->neq;
	nel = Parameters->nel;
	M2_out = FemStructs->AuxBuild->M2;
	R2_out = FemStructs->AuxBuild->R2;
	R1_out = FemStructs->AuxBuild->R1;
	invN2_out = FemStructs->AuxBuild->invN2;
	uB = FemStructs->uB;
	duB = FemStructs->duB;
	
	dzero(neq+1, F); // F do sistema M*Da = F*
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
		// Density
		//double theta[3];
		//FemFunctions->BC_General_theta(J1,J2,J3,Node,theta,FemFunctions->BC_theta);

		//Depois desse ponto vc utilizar os valores de theta[0],  theta[1] ou theta[2] convenientimente

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

		/*---------------------------------------------ALTERACAO Y^-1)-------------------------------------------------------*/
		double A0[4][4];		
		delta = FemFunctions->ShockCapture(Parameters, delta_old, gradUx, gradUy, Ax, Ay, A0, dUb, y23, y31, y12, x32, x13, x21, twoArea, e, Parameters->invY, Ub);

		//delta_soma += delta;
		//if(delta > delta_max)
		//	delta_max = delta;
		
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
		
		// *** Matriz de Massa MhB 12x4 e MBh 4x12, onde os coeficientes são iguais, assim calculamos apenas MhB
		double MhB;
		MhB = (3.0 * Area) / 20.0;
		
		// *** Matriz de Massa MBB 4x4
		double MBB;
		MBB = (81.0 * Area) / 280.0;
		
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
		

		deltaNMV = FemFunctions->ShockCapture_Micro(Parameters, deltaNMV_old, gradUx, gradUy, Ax, Ay, A0, dUb, y23, y31, y12, x32, x13, x21, twoArea, e, Parameters->invY, Ub);
	//	deltaNMV = delta;


		// *** Matriz de Rigidez KBB 4x4
		double KBB;
		
		//delta = delta/Cstab;
		KBB = (81.0*deltaNMV*(y23*y23 + y31*y31 + y12*y12 + y23*y31 + y23*y12 + y31*y12 + x32*x32 + x13*x13 + x21*x21 + x32*x13 + x32*x21 + x13*x21))/(40.0*Area);
		
		// Calculo auxiliar
		double alpha_dt = alpha*delta_t;
		
		// M1 = Mhh + alpha*dt*Khh
		double M1[12][12];
		// no1 - no1 
		M1[0][0] = Mhh1 + alpha_dt*(Khhg11 + Khhdd11);
		M1[0][1] =      + alpha_dt*Khhg12;
		M1[0][2] =      + alpha_dt*Khhg13;
		M1[0][3] =      + alpha_dt*Khhg14;
		M1[1][0] =      + alpha_dt*Khhg21;
		M1[1][1] = Mhh1 + alpha_dt*(Khhg22 + Khhdd11);
		M1[1][2] =      + alpha_dt*Khhg23;
		M1[1][3] =      + alpha_dt*Khhg24;
		M1[2][0] =      + alpha_dt*Khhg31;
		M1[2][1] =      + alpha_dt*Khhg32;
		M1[2][2] = Mhh1 + alpha_dt*(Khhg33 + Khhdd11);
		M1[2][3] =      + alpha_dt*Khhg34;
		M1[3][0] =      + alpha_dt*Khhg41;
		M1[3][1] =      + alpha_dt*Khhg42;
		M1[3][2] =      + alpha_dt*Khhg43;
		M1[3][3] = Mhh1 + alpha_dt*(Khhg44 + Khhdd11);
		
		// no1 - no2 
		M1[0][4] = Mhh2 + alpha_dt*(Khhg15 + Khhdd15);
		M1[0][5] =      + alpha_dt*Khhg16;
		M1[0][6] =      + alpha_dt*Khhg17;
		M1[0][7] =      + alpha_dt*Khhg18;
		M1[1][4] =      + alpha_dt*Khhg25;
		M1[1][5] = Mhh2 + alpha_dt*(Khhg26 + Khhdd15);
		M1[1][6] =      + alpha_dt*Khhg27;
		M1[1][7] =      + alpha_dt*Khhg28;
		M1[2][4] =      + alpha_dt*Khhg35;
		M1[2][5] =      + alpha_dt*Khhg36;
		M1[2][6] = Mhh2 + alpha_dt*(Khhg37 + Khhdd15);
		M1[2][7] =      + alpha_dt*Khhg38;
		M1[3][4] =      + alpha_dt*Khhg45;
		M1[3][5] =      + alpha_dt*Khhg46;
		M1[3][6] =      + alpha_dt*Khhg47;
		M1[3][7] = Mhh2 + alpha_dt*(Khhg48 + Khhdd15);
		
		// no1 - no3 
		M1[0][8]  = Mhh2 + alpha_dt*(Khhg19 + Khhdd19); 
		M1[0][9]  =      + alpha_dt*Khhg110;
		M1[0][10] =      + alpha_dt*Khhg111;
		M1[0][11] =      + alpha_dt*Khhg112;
		M1[1][8]  =      + alpha_dt*Khhg29;
		M1[1][9]  = Mhh2 + alpha_dt*(Khhg210 + Khhdd19);
		M1[1][10] =      + alpha_dt*Khhg211;
		M1[1][11] =      + alpha_dt*Khhg212;
		M1[2][8]  =      + alpha_dt*Khhg39;
		M1[2][9]  =      + alpha_dt*Khhg310;
		M1[2][10] = Mhh2 + alpha_dt*(Khhg311 + Khhdd19);
		M1[2][11] =      + alpha_dt*Khhg312;
		M1[3][8]  =      + alpha_dt*Khhg49;
		M1[3][9]  =      + alpha_dt*Khhg410;
		M1[3][10] =      + alpha_dt*Khhg411;
		M1[3][11] = Mhh2 + alpha_dt*(Khhg412 + Khhdd19);
		
		// no2 - no1 
		M1[4][0] = Mhh2 + alpha_dt*(Khhg11 + Khhdd15);
		M1[4][1] =      + alpha_dt*Khhg12;
		M1[4][2] =      + alpha_dt*Khhg13;
		M1[4][3] =      + alpha_dt*Khhg14;
		M1[5][0] =      + alpha_dt*Khhg21;
		M1[5][1] = Mhh2 + alpha_dt*(Khhg22 + Khhdd15);
		M1[5][2] =      + alpha_dt*Khhg23;
		M1[5][3] =      + alpha_dt*Khhg24;
		M1[6][0] =      + alpha_dt*Khhg31;
		M1[6][1] =      + alpha_dt*Khhg32;
		M1[6][2] = Mhh2 + alpha_dt*(Khhg33 + Khhdd15);
		M1[6][3] =      + alpha_dt*Khhg34;
		M1[7][0] =      + alpha_dt*Khhg41;
		M1[7][1] =      + alpha_dt*Khhg42;
		M1[7][2] =      + alpha_dt*Khhg43;
		M1[7][3] = Mhh2 + alpha_dt*(Khhg44 + Khhdd15);
		
		// no2 - no2 
		M1[4][4] = Mhh1 + alpha_dt*(Khhg15 + Khhdd55);
		M1[4][5] =      + alpha_dt*Khhg16;
		M1[4][6] =      + alpha_dt*Khhg17;
		M1[4][7] =      + alpha_dt*Khhg18;
		M1[5][4] =      + alpha_dt*Khhg25;
		M1[5][5] = Mhh1 + alpha_dt*(Khhg26 + Khhdd55);
		M1[5][6] =      + alpha_dt*Khhg27;
		M1[5][7] =      + alpha_dt*Khhg28;
		M1[6][4] =      + alpha_dt*Khhg35;
		M1[6][5] =      + alpha_dt*Khhg36;
		M1[6][6] = Mhh1 + alpha_dt*(Khhg37 + Khhdd55);
		M1[6][7] =      + alpha_dt*Khhg38;
		M1[7][4] =      + alpha_dt*Khhg45;
		M1[7][5] =      + alpha_dt*Khhg46;
		M1[7][6] =      + alpha_dt*Khhg47;
		M1[7][7] = Mhh1 + alpha_dt*(Khhg48 + Khhdd55);
		
		// no2 - no3 
		M1[4][8]  = Mhh2 + alpha_dt*(Khhg19 + Khhdd59); 
		M1[4][9]  =      + alpha_dt*Khhg110;
		M1[4][10] =      + alpha_dt*Khhg111;
		M1[4][11] =      + alpha_dt*Khhg112;
		M1[5][8]  =      + alpha_dt*Khhg29;
		M1[5][9]  = Mhh2 + alpha_dt*(Khhg210 + Khhdd59);
		M1[5][10] =      + alpha_dt*Khhg211;
		M1[5][11] =      + alpha_dt*Khhg212;
		M1[6][8]  =      + alpha_dt*Khhg39;
		M1[6][9]  =      + alpha_dt*Khhg310;
		M1[6][10] = Mhh2 + alpha_dt*(Khhg311 + Khhdd59);
		M1[6][11] =      + alpha_dt*Khhg312;
		M1[7][8]  =      + alpha_dt*Khhg49;
		M1[7][9]  =      + alpha_dt*Khhg410;
		M1[7][10] =      + alpha_dt*Khhg411;
		M1[7][11] = Mhh2 + alpha_dt*(Khhg412 + Khhdd59);
		
		// no3 - no1 
		M1[8][0]  = Mhh2 + alpha_dt*(Khhg11 + Khhdd19);
		M1[8][1]  =      + alpha_dt*Khhg12;
		M1[8][2]  =      + alpha_dt*Khhg13;
		M1[8][3]  =      + alpha_dt*Khhg14;
		M1[9][0]  =      + alpha_dt*Khhg21;
		M1[9][1]  = Mhh2 + alpha_dt*(Khhg22 + Khhdd19);
		M1[9][2]  =      + alpha_dt*Khhg23;
		M1[9][3]  =      + alpha_dt*Khhg24;
		M1[10][0] =      + alpha_dt*Khhg31;
		M1[10][1] =      + alpha_dt*Khhg32;
		M1[10][2] = Mhh2 + alpha_dt*(Khhg33 + Khhdd19);
		M1[10][3] =      + alpha_dt*Khhg34;
		M1[11][0] =      + alpha_dt*Khhg41;
		M1[11][1] =      + alpha_dt*Khhg42;
		M1[11][2] =      + alpha_dt*Khhg43;
		M1[11][3] = Mhh2 + alpha_dt*(Khhg44 + Khhdd19);
		
		// no3 - no2 
		M1[8][4]  = Mhh2 + alpha_dt*(Khhg15 + Khhdd59);
		M1[8][5]  =      + alpha_dt*Khhg16;
		M1[8][6]  =      + alpha_dt*Khhg17;
		M1[8][7]  =      + alpha_dt*Khhg18;
		M1[9][4]  =      + alpha_dt*Khhg25;
		M1[9][5]  = Mhh2 + alpha_dt*(Khhg26 + Khhdd59);
		M1[9][6]  =      + alpha_dt*Khhg27;
		M1[9][7]  =      + alpha_dt*Khhg28;
		M1[10][4] =      + alpha_dt*Khhg35;
		M1[10][5] =      + alpha_dt*Khhg36;
		M1[10][6] = Mhh2 + alpha_dt*(Khhg37 + Khhdd59);
		M1[10][7] =      + alpha_dt*Khhg38;
		M1[11][4] =      + alpha_dt*Khhg45;
		M1[11][5] =      + alpha_dt*Khhg46;
		M1[11][6] =      + alpha_dt*Khhg47;
		M1[11][7] = Mhh2 + alpha_dt*(Khhg48 + Khhdd59);
		
		// no3 - no3 
		M1[8][8]   = Mhh1 + alpha_dt*(Khhg19 + Khhdd99); 
		M1[8][9]   =      + alpha_dt*Khhg110;
		M1[8][10]  =      + alpha_dt*Khhg111;
		M1[8][11]  =      + alpha_dt*Khhg112;
		M1[9][8]   =      + alpha_dt*Khhg29;
		M1[9][9]   = Mhh1 + alpha_dt*(Khhg210 + Khhdd99);
		M1[9][10]  =      + alpha_dt*Khhg211;
		M1[9][11]  =      + alpha_dt*Khhg212;
		M1[10][8]  =      + alpha_dt*Khhg39; 
		M1[10][9]  =      + alpha_dt*Khhg310;
		M1[10][10] = Mhh1 + alpha_dt*(Khhg311 + Khhdd99);
		M1[10][11] =      + alpha_dt*Khhg312;
		M1[11][8]  =      + alpha_dt*Khhg49;
		M1[11][9]  =      + alpha_dt*Khhg410;
		M1[11][10] =      + alpha_dt*Khhg411;
		M1[11][11] = Mhh1 + alpha_dt*(Khhg412 + Khhdd99);
		
		// N1 = MhB + alpha*dt*KhB
		double N1[12][4];
		
		N1[0][0] = MhB + alpha_dt*KhB11;
		N1[0][1] =       alpha_dt*KhB12;
		N1[0][2] =       alpha_dt*KhB13;
		N1[0][3] =       alpha_dt*KhB14;
		N1[1][0] =       alpha_dt*KhB21;
		N1[1][1] = MhB + alpha_dt*KhB22;
		N1[1][2] =       alpha_dt*KhB23;
		N1[1][3] =       alpha_dt*KhB24; 
		N1[2][0] =       alpha_dt*KhB31;
		N1[2][1] =       alpha_dt*KhB32;
		N1[2][2] = MhB + alpha_dt*KhB33;
		N1[2][3] =       alpha_dt*KhB34;
		N1[3][0] =       alpha_dt*KhB41;
		N1[3][1] =       alpha_dt*KhB42;
		N1[3][2] =       alpha_dt*KhB43;
		N1[3][3] = MhB + alpha_dt*KhB44;
		
		N1[4][0] = MhB + alpha_dt*KhB51;
		N1[4][1] =       alpha_dt*KhB52;
		N1[4][2] =       alpha_dt*KhB53;
		N1[4][3] =       alpha_dt*KhB54;
		N1[5][0] =       alpha_dt*KhB61;
		N1[5][1] = MhB + alpha_dt*KhB62;
		N1[5][2] =       alpha_dt*KhB63;
		N1[5][3] =       alpha_dt*KhB64;
		N1[6][0] =       alpha_dt*KhB71;
		N1[6][1] =       alpha_dt*KhB72;
		N1[6][2] = MhB + alpha_dt*KhB73;
		N1[6][3] =       alpha_dt*KhB74;
		N1[7][0] =       alpha_dt*KhB81;
		N1[7][1] =       alpha_dt*KhB82;
		N1[7][2] =       alpha_dt*KhB83;
		N1[7][3] = MhB + alpha_dt*KhB84;
		
		N1[8][0]  = MhB + alpha_dt*KhB91;
		N1[8][1]  =       alpha_dt*KhB92;
		N1[8][2]  =       alpha_dt*KhB93;
		N1[8][3]  =       alpha_dt*KhB94;
		N1[9][0]  =       alpha_dt*KhB101;
		N1[9][1]  = MhB + alpha_dt*KhB102;
		N1[9][2]  =       alpha_dt*KhB103;
		N1[9][3]  =       alpha_dt*KhB104;
		N1[10][0] =       alpha_dt*KhB111;
		N1[10][1] =       alpha_dt*KhB112;
		N1[10][2] = MhB + alpha_dt*KhB113;
		N1[10][3] =       alpha_dt*KhB114;
		N1[11][0] =       alpha_dt*KhB121;
		N1[11][1] =       alpha_dt*KhB122;
		N1[11][2] =       alpha_dt*KhB123;
		N1[11][3] = MhB + alpha_dt*KhB124;
		
		//M2 = MBh + alpha*dt*KBh
		double M2[4][12];
		
		M2[0][0] = MhB + alpha_dt*KBh11;
		M2[0][1] =       alpha_dt*KBh12;
		M2[0][2] =       alpha_dt*KBh13;
		M2[0][3] =       alpha_dt*KBh14;
		M2[1][0] =       alpha_dt*KBh21;
		M2[1][1] = MhB + alpha_dt*KBh22;
		M2[1][2] =       alpha_dt*KBh23;
		M2[1][3] =       alpha_dt*KBh24; 
		M2[2][0] =       alpha_dt*KBh31;
		M2[2][1] =       alpha_dt*KBh32;
		M2[2][2] = MhB + alpha_dt*KBh33;
		M2[2][3] =       alpha_dt*KBh34;
		M2[3][0] =       alpha_dt*KBh41;
		M2[3][1] =       alpha_dt*KBh42;
		M2[3][2] =       alpha_dt*KBh43;
		M2[3][3] = MhB + alpha_dt*KBh44;
		
		M2[0][4] = MhB + alpha_dt*KBh15;
		M2[0][5] =       alpha_dt*KBh16;
		M2[0][6] =       alpha_dt*KBh17;
		M2[0][7] =       alpha_dt*KBh18;
		M2[1][4] =       alpha_dt*KBh25;
		M2[1][5] = MhB + alpha_dt*KBh26;
		M2[1][6] =       alpha_dt*KBh27;
		M2[1][7] =       alpha_dt*KBh28;
		M2[2][4] =       alpha_dt*KBh35;
		M2[2][5] =       alpha_dt*KBh36;
		M2[2][6] = MhB + alpha_dt*KBh37;
		M2[2][7] =       alpha_dt*KBh38;
		M2[3][4] =       alpha_dt*KBh45;
		M2[3][5] =       alpha_dt*KBh46;
		M2[3][6] =       alpha_dt*KBh47;
		M2[3][7] = MhB + alpha_dt*KBh48;
		
		M2[0][8]  = MhB + alpha_dt*KBh19;
		M2[0][9]  =       alpha_dt*KBh110;
		M2[0][10] =       alpha_dt*KBh111;
		M2[0][11] =       alpha_dt*KBh112;
		M2[1][8]  =       alpha_dt*KBh29;
		M2[1][9]  = MhB + alpha_dt*KBh210;
		M2[1][10] =       alpha_dt*KBh211;
		M2[1][11] =       alpha_dt*KBh212;
		M2[2][8]  =       alpha_dt*KBh39;
		M2[2][9]  =       alpha_dt*KBh310;
		M2[2][10] = MhB + alpha_dt*KBh311;
		M2[2][11] =       alpha_dt*KBh312;
		M2[3][8]  =       alpha_dt*KBh49;
		M2[3][9]  =       alpha_dt*KBh410;
		M2[3][10] =       alpha_dt*KBh411;
		M2[3][11] = MhB + alpha_dt*KBh412;
		
		// N2^{-1} = (MBB + alpha*dt*KBB)^{-1}
		double invN2 = 1.0 / (MBB + alpha_dt*KBB);
		/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mhh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
		double Mhh[12][12];
		// no1 - no1 
		Mhh[0][0] = Mhh1;
		Mhh[0][1] = 0.0;
		Mhh[0][2] = 0.0;
		Mhh[0][3] = 0.0;
		Mhh[1][0] = 0.0;
		Mhh[1][1] = Mhh1;
		Mhh[1][2] = 0.0;
		Mhh[1][3] = 0.0;
		Mhh[2][0] = 0.0;
		Mhh[2][1] = 0.0;
		Mhh[2][2] = Mhh1;
		Mhh[2][3] = 0.0;
		Mhh[3][0] = 0.0;
		Mhh[3][1] = 0.0;
		Mhh[3][2] = 0.0;
		Mhh[3][3] = Mhh1;
		
		// no1 - no2 
		Mhh[0][4] = Mhh2;
		Mhh[0][5] = 0.0;
		Mhh[0][6] = 0.0;
		Mhh[0][7] = 0.0;
		Mhh[1][4] = 0.0;
		Mhh[1][5] = Mhh2;
		Mhh[1][6] = 0.0;
		Mhh[1][7] = 0.0;
		Mhh[2][4] = 0.0;
		Mhh[2][5] = 0.0;
		Mhh[2][6] = Mhh2;
		Mhh[2][7] = 0.0;
		Mhh[3][4] = 0.0;
		Mhh[3][5] = 0.0;
		Mhh[3][6] = 0.0;
		Mhh[3][7] = Mhh2;
		
		// no1 - no3 
		Mhh[0][8]  = Mhh2; 
		Mhh[0][9]  = 0.0;
		Mhh[0][10] = 0.0;
		Mhh[0][11] = 0.0;
		Mhh[1][8]  = 0.0;
		Mhh[1][9]  = Mhh2;
		Mhh[1][10] = 0.0;
		Mhh[1][11] = 0.0;
		Mhh[2][8]  = 0.0;
		Mhh[2][9]  = 0.0;
		Mhh[2][10] = Mhh2;
		Mhh[2][11] = 0.0;
		Mhh[3][8]  = 0.0;
		Mhh[3][9]  = 0.0;
		Mhh[3][10] = 0.0;
		Mhh[3][11] = Mhh2;
		
		// no2 - no1 
		Mhh[4][0] = Mhh2;
		Mhh[4][1] = 0.0;
		Mhh[4][2] = 0.0;
		Mhh[4][3] = 0.0;
		Mhh[5][0] = 0.0;
		Mhh[5][1] = Mhh2;
		Mhh[5][2] = 0.0;
		Mhh[5][3] = 0.0;
		Mhh[6][0] = 0.0;
		Mhh[6][1] = 0.0;
		Mhh[6][2] = Mhh2;
		Mhh[6][3] = 0.0;
		Mhh[7][0] = 0.0;
		Mhh[7][1] = 0.0;
		Mhh[7][2] = 0.0;
		Mhh[7][3] = Mhh2;
		
		// no2 - no2 
		Mhh[4][4] = Mhh1;
		Mhh[4][5] = 0.0;
		Mhh[4][6] = 0.0;
		Mhh[4][7] = 0.0;
		Mhh[5][4] = 0.0;
		Mhh[5][5] = Mhh1;
		Mhh[5][6] = 0.0;
		Mhh[5][7] = 0.0;
		Mhh[6][4] = 0.0;
		Mhh[6][5] = 0.0;
		Mhh[6][6] = Mhh1;
		Mhh[6][7] = 0.0;
		Mhh[7][4] = 0.0;
		Mhh[7][5] = 0.0;
		Mhh[7][6] = 0.0;
		Mhh[7][7] = Mhh1;
		
		// no2 - no3 
		Mhh[4][8]  = Mhh2; 
		Mhh[4][9]  = 0.0;
		Mhh[4][10] = 0.0;
		Mhh[4][11] = 0.0;
		Mhh[5][8]  = 0.0;
		Mhh[5][9]  = Mhh2;
		Mhh[5][10] = 0.0;
		Mhh[5][11] = 0.0;
		Mhh[6][8]  = 0.0;
		Mhh[6][9]  = 0.0;
		Mhh[6][10] = Mhh2;
		Mhh[6][11] = 0.0;
		Mhh[7][8]  = 0.0;
		Mhh[7][9]  = 0.0;
		Mhh[7][10] = 0.0;
		Mhh[7][11] = Mhh2;
		
		// no3 - no1 
		Mhh[8][0]  = Mhh2;
		Mhh[8][1]  = 0.0;
		Mhh[8][2]  = 0.0;
		Mhh[8][3]  = 0.0;
		Mhh[9][0]  = 0.0;
		Mhh[9][1]  = Mhh2;
		Mhh[9][2]  = 0.0;
		Mhh[9][3]  = 0.0;
		Mhh[10][0] = 0.0;
		Mhh[10][1] = 0.0;
		Mhh[10][2] = Mhh2;
		Mhh[10][3] = 0.0;
		Mhh[11][0] = 0.0;
		Mhh[11][1] = 0.0;
		Mhh[11][2] = 0.0;
		Mhh[11][3] = Mhh2;
		
		// no3 - no2 
		Mhh[8][4]  = Mhh2;
		Mhh[8][5]  = 0.0;
		Mhh[8][6]  = 0.0;
		Mhh[8][7]  = 0.0;
		Mhh[9][4]  = 0.0;
		Mhh[9][5]  = Mhh2;
		Mhh[9][6]  = 0.0;
		Mhh[9][7]  = 0.0;
		Mhh[10][4] = 0.0;
		Mhh[10][5] = 0.0;
		Mhh[10][6] = Mhh2;
		Mhh[10][7] = 0.0;
		Mhh[11][4] = 0.0;
		Mhh[11][5] = 0.0;
		Mhh[11][6] = 0.0;
		Mhh[11][7] = Mhh2;
		
		// no3 - no3 
		Mhh[8][8]   = Mhh1; 
		Mhh[8][9]   = 0.0;
		Mhh[8][10]  = 0.0;
		Mhh[8][11]  = 0.0;
		Mhh[9][8]   = 0.0;
		Mhh[9][9]   = Mhh1;
		Mhh[9][10]  = 0.0;
		Mhh[9][11]  = 0.0;
		Mhh[10][8]  = 0.0; 
		Mhh[10][9]  = 0.0;
		Mhh[10][10] = Mhh1;
		Mhh[10][11] = 0.0;
		Mhh[11][8]  = 0.0;
		Mhh[11][9]  = 0.0;
		Mhh[11][10] = 0.0;
		Mhh[11][11] = Mhh1;
		/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Khh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
		double Khh[12][12];
		// no1 - no1 
		Khh[0][0] = (Khhg11 + Khhdd11);
		Khh[0][1] = Khhg12;
		Khh[0][2] = Khhg13;
		Khh[0][3] = Khhg14;
		Khh[1][0] = Khhg21;
		Khh[1][1] = (Khhg22 + Khhdd11);
		Khh[1][2] = Khhg23;
		Khh[1][3] = Khhg24;
		Khh[2][0] = Khhg31;
		Khh[2][1] = Khhg32;
		Khh[2][2] = (Khhg33 + Khhdd11);
		Khh[2][3] = Khhg34;
		Khh[3][0] = Khhg41;
		Khh[3][1] = Khhg42;
		Khh[3][2] = Khhg43;
		Khh[3][3] = (Khhg44 + Khhdd11);
		
		// no1 - no2 
		Khh[0][4] = (Khhg15 + Khhdd15);
		Khh[0][5] = Khhg16;
		Khh[0][6] = Khhg17;
		Khh[0][7] = Khhg18;
		Khh[1][4] = Khhg25;
		Khh[1][5] = (Khhg26 + Khhdd15);
		Khh[1][6] = Khhg27;
		Khh[1][7] = Khhg28;
		Khh[2][4] = Khhg35;
		Khh[2][5] = Khhg36;
		Khh[2][6] = (Khhg37 + Khhdd15);
		Khh[2][7] = Khhg38;
		Khh[3][4] = Khhg45;
		Khh[3][5] = Khhg46;
		Khh[3][6] = Khhg47;
		Khh[3][7] = (Khhg48 + Khhdd15);
		
		// no1 - no3 
		Khh[0][8]  = (Khhg19 + Khhdd19); 
		Khh[0][9]  = Khhg110;
		Khh[0][10] = Khhg111;
		Khh[0][11] = Khhg112;
		Khh[1][8]  = Khhg29;
		Khh[1][9]  = (Khhg210 + Khhdd19);
		Khh[1][10] = Khhg211;
		Khh[1][11] = Khhg212;
		Khh[2][8]  = Khhg39;
		Khh[2][9]  = Khhg310;
		Khh[2][10] = (Khhg311 + Khhdd19);
		Khh[2][11] = Khhg312;
		Khh[3][8]  = Khhg49;
		Khh[3][9]  = Khhg410;
		Khh[3][10] = Khhg411;
		Khh[3][11] = (Khhg412 + Khhdd19);
		
		// no2 - no1 
		Khh[4][0] = (Khhg11 + Khhdd15);
		Khh[4][1] = Khhg12;
		Khh[4][2] = Khhg13;
		Khh[4][3] = Khhg14;
		Khh[5][0] = Khhg21;
		Khh[5][1] = (Khhg22 + Khhdd15);
		Khh[5][2] = Khhg23;
		Khh[5][3] = Khhg24;
		Khh[6][0] = Khhg31;
		Khh[6][1] = Khhg32;
		Khh[6][2] = (Khhg33 + Khhdd15);
		Khh[6][3] = Khhg34;
		Khh[7][0] = Khhg41;
		Khh[7][1] = Khhg42;
		Khh[7][2] = Khhg43;
		Khh[7][3] = (Khhg44 + Khhdd15);
		
		// no2 - no2 
		Khh[4][4] = (Khhg15 + Khhdd55);
		Khh[4][5] = Khhg16;
		Khh[4][6] = Khhg17;
		Khh[4][7] = Khhg18;
		Khh[5][4] = Khhg25;
		Khh[5][5] = (Khhg26 + Khhdd55);
		Khh[5][6] = Khhg27;
		Khh[5][7] = Khhg28;
		Khh[6][4] = Khhg35;
		Khh[6][5] = Khhg36;
		Khh[6][6] = (Khhg37 + Khhdd55);
		Khh[6][7] = Khhg38;
		Khh[7][4] = Khhg45;
		Khh[7][5] = Khhg46;
		Khh[7][6] = Khhg47;
		Khh[7][7] = (Khhg48 + Khhdd55);
		
		// no2 - no3 
		Khh[4][8]  = (Khhg19 + Khhdd59); 
		Khh[4][9]  = Khhg110;
		Khh[4][10] = Khhg111;
		Khh[4][11] = Khhg112;
		Khh[5][8]  = Khhg29;
		Khh[5][9]  = (Khhg210 + Khhdd59);
		Khh[5][10] = Khhg211;
		Khh[5][11] = Khhg212;
		Khh[6][8]  = Khhg39;
		Khh[6][9]  = Khhg310;
		Khh[6][10] = (Khhg311 + Khhdd59);
		Khh[6][11] = Khhg312;
		Khh[7][8]  = Khhg49;
		Khh[7][9]  = Khhg410;
		Khh[7][10] = Khhg411;
		Khh[7][11] = (Khhg412 + Khhdd59);
		
		// no3 - no1 
		Khh[8][0]  = (Khhg11 + Khhdd19);
		Khh[8][1]  = Khhg12;
		Khh[8][2]  = Khhg13;
		Khh[8][3]  = Khhg14;
		Khh[9][0]  = Khhg21;
		Khh[9][1]  = (Khhg22 + Khhdd19);
		Khh[9][2]  = Khhg23;
		Khh[9][3]  = Khhg24;
		Khh[10][0] = Khhg31;
		Khh[10][1] = Khhg32;
		Khh[10][2] = (Khhg33 + Khhdd19);
		Khh[10][3] = Khhg34;
		Khh[11][0] = Khhg41;
		Khh[11][1] = Khhg42;
		Khh[11][2] = Khhg43;
		Khh[11][3] = (Khhg44 + Khhdd19);
		
		// no3 - no2 
		Khh[8][4]  = (Khhg15 + Khhdd59);
		Khh[8][5]  = Khhg16;
		Khh[8][6]  = Khhg17;
		Khh[8][7]  = Khhg18;
		Khh[9][4]  = Khhg25;
		Khh[9][5]  = (Khhg26 + Khhdd59);
		Khh[9][6]  = Khhg27;
		Khh[9][7]  = Khhg28;
		Khh[10][4] = Khhg35;
		Khh[10][5] = Khhg36;
		Khh[10][6] = (Khhg37 + Khhdd59);
		Khh[10][7] = Khhg38;
		Khh[11][4] = Khhg45;
		Khh[11][5] = Khhg46;
		Khh[11][6] = Khhg47;
		Khh[11][7] = (Khhg48 + Khhdd59);
		
		// no3 - no3 
		Khh[8][8]   = (Khhg19 + Khhdd99); 
		Khh[8][9]   = Khhg110;
		Khh[8][10]  = Khhg111;
		Khh[8][11]  = Khhg112;
		Khh[9][8]   = Khhg29;
		Khh[9][9]   = (Khhg210 + Khhdd99);
		Khh[9][10]  = Khhg211;
		Khh[9][11]  = Khhg212;
		Khh[10][8]  = Khhg39; 
		Khh[10][9]  = Khhg310;
		Khh[10][10] = (Khhg311 + Khhdd99);
		Khh[10][11] = Khhg312;
		Khh[11][8]  = Khhg49;
		Khh[11][9]  = Khhg410;
		Khh[11][10] = Khhg411;
		Khh[11][11] = (Khhg412 + Khhdd99);

		/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MBh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
		double MBh[4][12];
		
		MBh[0][0] = MhB;
		MBh[0][1] = 0.0;
		MBh[0][2] = 0.0;
		MBh[0][3] = 0.0;
		MBh[1][0] = 0.0;
		MBh[1][1] = MhB;
		MBh[1][2] = 0.0;
		MBh[1][3] = 0.0;
		MBh[2][0] = 0.0;
		MBh[2][1] = 0.0;
		MBh[2][2] = MhB;
		MBh[2][3] = 0.0;
		MBh[3][0] = 0.0;
		MBh[3][1] = 0.0;
		MBh[3][2] = 0.0;
		MBh[3][3] = MhB;
		
		MBh[0][4] = MhB;
		MBh[0][5] = 0.0;
		MBh[0][6] = 0.0;
		MBh[0][7] = 0.0;
		MBh[1][4] = 0.0;
		MBh[1][5] = MhB;
		MBh[1][6] = 0.0;
		MBh[1][7] = 0.0;
		MBh[2][4] = 0.0;
		MBh[2][5] = 0.0;
		MBh[2][6] = MhB;
		MBh[2][7] = 0.0;
		MBh[3][4] = 0.0;
		MBh[3][5] = 0.0;
		MBh[3][6] = 0.0;
		MBh[3][7] = MhB;
		
		MBh[0][8]  = MhB;
		MBh[0][9]  = 0.0;
		MBh[0][10] = 0.0;
		MBh[0][11] = 0.0;
		MBh[1][8]  = 0.0;
		MBh[1][9]  = MhB;
		MBh[1][10] = 0.0;
		MBh[1][11] = 0.0;
		MBh[2][8]  = 0.0;
		MBh[2][9]  = 0.0;
		MBh[2][10] = MhB;
		MBh[2][11] = 0.0;
		MBh[3][8]  = 0.0;
		MBh[3][9]  = 0.0;
		MBh[3][10] = 0.0;
		MBh[3][11] = MhB;
		/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KBh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
		double KBh[4][12];
	
		KBh[0][0] = KBh11; 
		KBh[0][1] = KBh12;
		KBh[0][2] = KBh13;
		KBh[0][3] = KBh14;
		KBh[1][0] = KBh21; 
		KBh[1][1] = KBh22;
		KBh[1][2] = KBh23;
		KBh[1][3] = KBh24;
		KBh[2][0] = KBh31; 
		KBh[2][1] = KBh32;
		KBh[2][2] = KBh33;
		KBh[2][3] = KBh34;
		KBh[3][0] = KBh41; 
		KBh[3][1] = KBh42;
		KBh[3][2] = KBh43;
		KBh[3][3] = KBh44;

		KBh[0][4] = KBh15; 
		KBh[0][5] = KBh16;
		KBh[0][6] = KBh17;
		KBh[0][7] = KBh18;
		KBh[1][4] = KBh25; 
		KBh[1][5] = KBh26;
		KBh[1][6] = KBh27;
		KBh[1][7] = KBh28;
		KBh[2][4] = KBh35; 
		KBh[2][5] = KBh36;
		KBh[2][6] = KBh37;
		KBh[2][7] = KBh38;
		KBh[3][4] = KBh45; 
		KBh[3][5] = KBh46;
		KBh[3][6] = KBh47;
		KBh[3][7] = KBh48;

		KBh[0][8] = KBh19; 
		KBh[0][9] = KBh110;
		KBh[0][10] = KBh111;
		KBh[0][11] = KBh112;
		KBh[1][8] = KBh29; 
		KBh[1][9] = KBh210;
		KBh[1][10] = KBh211;
		KBh[1][11] = KBh212;
		KBh[2][8] = KBh39; 
		KBh[2][9] = KBh310;
		KBh[2][10] = KBh311;
		KBh[2][11] = KBh312;
		KBh[3][8] = KBh49; 
		KBh[3][9] = KBh410;
		KBh[3][10] = KBh411;
		KBh[3][11] = KBh412;
		
		/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KhB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
		double KhB[12][4];
		
		KhB[0][0] = KhB11;
		KhB[0][1] = KhB12;
		KhB[0][2] = KhB13;
		KhB[0][3] = KhB14;
		KhB[1][0] = KhB21;
		KhB[1][1] = KhB22;		
		KhB[1][2] = KhB23;
		KhB[1][3] = KhB24;
		KhB[2][0] = KhB31;
		KhB[2][1] = KhB32;
		KhB[2][2] = KhB33;
		KhB[2][3] = KhB34;
		KhB[3][0] = KhB41;
		KhB[3][1] = KhB42;
		KhB[3][2] = KhB43;
		KhB[3][3] = KhB44;
		
		KhB[4][0] = KhB51;
		KhB[4][1] = KhB52;
		KhB[4][2] = KhB53;
		KhB[4][3] = KhB54;
		KhB[5][0] = KhB61;
		KhB[5][1] = KhB62;		
		KhB[5][2] = KhB63;
		KhB[5][3] = KhB64;
		KhB[6][0] = KhB71;
		KhB[6][1] = KhB72;
		KhB[6][2] = KhB73;
		KhB[6][3] = KhB74;
		KhB[7][0] = KhB81;
		KhB[7][1] = KhB82;
		KhB[7][2] = KhB83;
		KhB[7][3] = KhB84;
						
		KhB[8][0] = KhB91;
		KhB[8][1] = KhB92;
		KhB[8][2] = KhB93;
		KhB[8][3] = KhB94;
		KhB[9][0] = KhB101;
		KhB[9][1] = KhB102;		
		KhB[9][2] = KhB103;
		KhB[9][3] = KhB104;
		KhB[10][0] = KhB111;
		KhB[10][1] = KhB112;
		KhB[10][2] = KhB113;
		KhB[10][3] = KhB114;
		KhB[11][0] = KhB121;
		KhB[11][1] = KhB122;
		KhB[11][2] = KhB123;
		KhB[11][3] = KhB124;		

		// Valor de uB e duB no elemento
		int eNDOF = e*NDOF;
		double uBaux[4], duBaux[4];
		uBaux[0] = uB[eNDOF];
		uBaux[1] = uB[eNDOF + 1];
		uBaux[2] = uB[eNDOF + 2];
		uBaux[3] = uB[eNDOF + 3];
		
		duBaux[0] = duB[eNDOF];
		duBaux[1] = duB[eNDOF + 1];
		duBaux[2] = duB[eNDOF + 2];
		duBaux[3] = duB[eNDOF + 3];
		
		//Montagem da matriz Me
		int size = NDOF*NNOEL;
		double soma;
		for(i = 0; i < size; i++){
			for(j = 0; j < size; j++){
				soma = 0.0;
				for (k = 0; k < NDOF; k++){
					soma += N1[i][k]*M2[k][j];
				}
				Me[i][j] = M1[i][j] - invN2*soma;
			}
		}
				
		double theta1, theta2, theta3; // Dessa forma não é possivel rotacionar os vetores R1 e R2
		if (Node[J1].v1Type == -1){
			theta1 = FemFunctions->BC_theta(x1,y1);
			rotationNMV(1, theta1, Me, Mhh, Khh, KhB, N1, MBh, KBh);
		}				
		if (Node[J2].v1Type == -1){
			theta2=FemFunctions->BC_theta(x2,y2);
			rotationNMV(2, theta2, Me, Mhh, Khh, KhB, N1, MBh, KBh);
		}
		if (Node[J3].v1Type == -1){
			theta3=FemFunctions->BC_theta(x3,y3);
			rotationNMV(3, theta3, Me, Mhh, Khh, KhB, N1, MBh, KBh);
		}						

		//Residuo
		double R1[12], R2[4], Res[12];
		// R1 = Fh - (Mhh*du + Khh*u) - (MhB*duB + KhB*uB)
		for (i=0; i<12; i++){
				R1[i] = 0.0;
				Res[i] = 0.0;				
				for (j=0; j<12; j++){
					R1[i] += -(Mhh[i][j]*dUe[j] + Khh[i][j]*Ue[j]);									
				}
				//Res[i] = R1[i];										
		}
		for (i=0; i<12; i++){								
			for (j=0; j<4; j++)
				R1[i] += -(MBh[j][i]*duBaux[j] + KhB[i][j]*uBaux[j]);  //(MBh)^T*dUB = MhB*dUB	
			Res[i] = R1[i];
		}
		//R2 = FB - (MBh*du + KBh*u) - (MBB*duB + KBB*uB)
		//R2 = - (MBh*du + KBh*u):
		for (i=0; i<4; i++){
			R2[i] = 0.0;
			for (j=0; j<12; j++)
				R2[i] += -(MBh[i][j]*dUe[j] + KBh[i][j]*Ue[j]);
		}

		//R2 = R2 - (MBB*duB + KBB*uB):
		R2[0]  += -(MBB*duBaux[0] + KBB*uBaux[0]);
		R2[1]  += -(MBB*duBaux[1] + KBB*uBaux[1]);
		R2[2]  += -(MBB*duBaux[2] + KBB*uBaux[2]);
		R2[3]  += -(MBB*duBaux[3] + KBB*uBaux[3]);

		// Montagem do vetor F 
		for (i=0; i<12; i++)	
			F[lm[e][i]] += R1[i]  - invN2*(N1[i][0]*R2[0]  + N1[i][1]*R2[1]  + N1[i][2]*R2[2]  + N1[i][3]*R2[3]);

		F[neq] = 0;


		// Matrizes e vetores para serem utilizados na funcao calcula_DaB		for (i=0; i<4; i++)
		for (i=0; i<4; i++)
			R2_out[e][i] = R2[i];

		for (i=0; i<12; i++)
			R1_out[e][i] = Res[i];

		invN2_out[e] = invN2;

		int k = 0;
		for (i=0; i<4; i++)
			for (j=0; j<12; j++, k++)
				M2_out[e][k] = M2[i][j];

		//Assembly of matrices Me
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, e, Me);

	}//for elemento
	//double	media_delta;
	//media_delta = delta_soma/nel;
	//printf("\n\n Media delta= %.12f\t delta maximo=%.12f\n\n",media_delta,delta_max);
		
	// Liberacao dos espacos alocados
	free(U);
	free(dU);
	

	return 0;
}// end build
















