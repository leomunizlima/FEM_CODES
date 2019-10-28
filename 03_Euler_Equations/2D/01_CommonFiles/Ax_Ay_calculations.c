#include "EulerEquations.h"

	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--- VLR ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void VLR_Preconditioner_Ax_Ay_calculations(double cv, double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4], double Dpmax)
{
	double rho = Ub[0];
	double v1 = Ub[1] / Ub[0];
	double v2 = Ub[2] / Ub[0];
	double E  = Ub[3] / Ub[0];
	double norma_U23 = (Ub[1] * Ub[1]) + (Ub[2] * Ub[2]);
	double norma_u2 = v1*v1 + v2*v2; 					//||u||^2
	double a1, a2, a3, a4, a5, a6, a7, a8, a9;
	double tau, beta, beta2, M, M2, Mbeta, R, c, c2;
	double P[4][4], AxNP[4][4], AyNP[4][4];				//AxNP e AyNP: Matrizes Ax e Ay antes do precondicionamento
	int i, j, k;

	c = sqrt(gamma*(gamma - 1.0)*(E - 0.5*norma_u2));	//velocidade do som no meio
	//if(norma_u2<1e-10)
	//	norma_u2 = 0.001;
	M = sqrt(norma_u2)/c;								//Mach local
	R = 0.4*cv;											//constante dos gases ideais

	//if(norma_u2<1e-6 || c < 1e-6)
	//	printf("\n\n||u||=%.12f\t c=%.12f\n\n",norma_u2,c);
	
	if(M >= 1.0 && M < 1.01)							//definicao de beta de acordo com Moragues (2016) com epsilon=0.01
		Mbeta = 1.01;
	else if(M > 0.99 && M < 1.0)
		Mbeta = 0.99;
	else
		Mbeta = M; 										//Não estou excluindo M=0! Mas não há problema.
	beta = sqrt(fabs(1.0 - Mbeta*Mbeta));
	
	if(M < 1.0)											//definicao de tau de acordo com Moragues (2016)
		tau = beta;
	else
		tau = beta/M;
	
	beta2 = beta*beta;
	
	M2 = M*M;
	c2 = c*c;
	//coeficientes para definicao da matriz P (VLR) - conferido!
	a1 = 1.0/norma_u2*(1.0 + tau/beta2 - tau) + 1.0/c2*( R/cv * (tau*(1.0 - M2)/beta2) - tau/beta2 + R/cv);		//a

	a2 = tau*(M2 - 1.0)/beta2*(1.0 + 0.5*(R/cv)*M2) - 0.5*(R/cv)*M2;						//b

	a3 = (R/cv)*(1.0/c2)*(tau*(M2 - 1.0)/beta2 - 1.0);								//c

	a4 = (R/cv)*(1.0/c2)*(1.0 - tau*M2/beta2) - tau/(c2*beta2);							//d

	a5 = 1.0 + tau*M2/beta2 + 0.5*(R/cv)*M2*(tau*M2/beta2 - 1.0);							//e

	a6 = (R/cv)*(1.0/c2)*(tau*M2/beta2 - 1.0);									//f
	
	a7 = 1.0 + (1.0 - cv/R)*tau/beta2 + (R/cv - 1.5)*tau*M2/beta2 + 0.5*(R/cv)*M2*(1.0 - tau*M2/beta2);		//g

	a8 = norma_u2*( (cv/R - 1.0 + M2*(1.0 - 0.5*R/cv))*tau/beta2 + 0.25*(R/cv)*M2*(tau*M2/beta2 - 1.0) - 0.5);	//h

	a9 = (1.0 - R/cv)*tau*M2/beta2 + 0.5*(R/cv)*M2*(tau*M2/beta2 - 1.0);						//i
	

	//Matriz P (VLR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	P[0][0] = a5;
	P[1][0] = a2*v1;	
	P[2][0] = a2*v2;	
	P[3][0] = a8;

	P[0][1] = a4*v1;
	P[1][1] = a1*v1*v1 + tau;		
	P[2][1] = a1*v1*v2;	
	P[3][1] = a7*v1;

	P[0][2] = a4*v2;
	P[1][2] = a1*v1*v2;	
	P[2][2] = a1*v2*v2 + tau;	
	P[3][2] = a7*v2;	

	P[0][3] = a6;
	P[1][3] = a3*v1;	
	P[2][3] = a3*v2;			
	P[3][3] = a9;

	// *** Ax coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	AxNP[0][0] = 0.0;
	AxNP[0][1] = 1.0;
	AxNP[0][2] = 0.0;
	AxNP[0][3] = 0.0;
	AxNP[1][0] = -(v1 * v1) + (((gamma - 1.0) * norma_U23 * 0.5) / (rho * rho));
	AxNP[1][1] = (3.0 - gamma) * v1;
	AxNP[1][2] = -(gamma - 1.0) * v2;
	AxNP[1][3] = gamma - 1.0;
	AxNP[2][0] = - v1*v2;
	AxNP[2][1] = v2;
	AxNP[2][2] = v1;
	AxNP[2][3] = 0.0;
	AxNP[3][0] = (-gamma * v1 * E) + ((gamma - 1.0) * v1 * (norma_U23 / (rho * rho)));
	AxNP[3][1] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v1 * v1)));
	AxNP[3][2] = -((gamma - 1.0) * v1 * v2);
	AxNP[3][3] = gamma * v1;

	// *** Ay coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	AyNP[0][0] = 0.0;
	AyNP[0][1] = 0.0;
	AyNP[0][2] = 1.0;
	AyNP[0][3] = 0.0;
	AyNP[1][0] = AxNP[2][0];
	AyNP[1][1] = v2;
	AyNP[1][2] = v1;
	AyNP[1][3] = 0.0;
	AyNP[2][0] = -(v2 * v2) + (((gamma - 1.0) * norma_U23 * 0.5) / (rho * rho)); 
	AyNP[2][1] = -(gamma - 1.0) * v1;
	AyNP[2][2] = (3.0 - gamma) * v2;
	AyNP[2][3] = AxNP[1][3];
	AyNP[3][0] = (- gamma * v2 * E) + (((gamma - 1.0) * v2 * norma_U23) / (rho * rho));
	AyNP[3][1] = AxNP[3][2];
	AyNP[3][2] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v2 * v2)));
	AyNP[3][3] = gamma * v2;
	
	//Produto PAxNP = Ax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for(i=0; i<4; i++)
		for(j=0; j<4; j++){
			Ax[i][j] = 0.0;
			for(k=0; k<4; k++)
				Ax[i][j] += P[i][k] * AxNP[k][j];
		}	
	//Produto PAy = Ay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	for(i=0; i<4; i++)
		for(j=0; j<4; j++){
			Ay[i][j] = 0.0;
			for(k=0; k<4; k++)
				Ay[i][j] += P[i][k] * AyNP[k][j];
		}	
}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--- WSCM ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void WSCM_Preconditioner_Ax_Ay_calculations(double cv, double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4], double Dpmax)
{
	double rho = Ub[0];
	double v1 = Ub[1] / Ub[0];
	double v2 = Ub[2] / Ub[0];
	double E  = Ub[3] / Ub[0];
	double norma_U23 = (Ub[1] * Ub[1]) + (Ub[2] * Ub[2]);
	double norma_u2 = v1*v1 + v2*v2; 							//||u||^2	
	double sigma, parameter_euler, epsilon, alpha, delta, theta, H, p, M2lim, M, M2, M2max, c, c2; 
	double P[4][4], AxNP[4][4], AyNP[4][4];						//AxNP e AyNP: Matrizes Ax e Ay antes do precondicionamento
	int i, j, k;
	//Parameters definitions according to Colin et al 2011
	p = (gamma - 1.0)*(rho*E - 0.5*rho*norma_u2);				//pressure
	c = sqrt(gamma*(gamma - 1.0)*(E - 0.5*norma_u2));			//velocidade do som no meio
	c2 = c*c;
	M = sqrt(norma_u2)/c;										//Mach local
	M2lim = 1e-5;
	M2 = M*M;
	delta = 0.5;											//free parameter: delta = 0 => WSCM=WS and delta = 1 => WSCM=CM
	sigma = 2.0;												
	parameter_euler =  sigma*Dpmax/(rho*c2);
	H = E + p/rho;												//Total enthalpy

	//%%%%%%%%%--epsilon definition--%%%%%%%%%%%%%%%%%%%%
	if(M2lim < M2)
		M2max = M2;
	else
		M2max = M2lim;
	if(M2max < parameter_euler){
		M2max = parameter_euler;
		//printf("\n\nParameter Euler:%lf", M2max);
		//getchar();
	}
	if(M2max<1.0)
		epsilon = M2max;
	else
		epsilon = 1.0;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	alpha = 1.0/c2*(gamma - 1.0)*( (1.0 - delta)*epsilon - 1.0 ); 
	
	theta = 0.5*norma_u2 + delta*(epsilon*c2)/((gamma - 1.0)*((1.0 - delta)*epsilon - 1.0)); 
			
	//Matriz P (WSCM) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	P[0][0] = alpha*theta + 1.0;
	P[1][0] = alpha*(v1*theta);	
	P[2][0] = alpha*(v2*theta);	
	P[3][0] = alpha*(H*theta);

	P[0][1] = -alpha*v1;
	P[1][1] = -alpha*v1*v1 + 1.0;		
	P[2][1] = -alpha*v1*v2;	
	P[3][1] = -alpha*v1*H;

	P[0][2] = -alpha*v2;
	P[1][2] = -alpha*v1*v2;	
	P[2][2] = -alpha*v2*v2 + 1.0;	
	P[3][2] = -alpha*v2*H;	

	P[0][3] = alpha;
	P[1][3] = alpha*v1;	
	P[2][3] = alpha*v2;			
	P[3][3] = alpha*H + 1.0;
		
	// *** Ax coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	AxNP[0][0] = 0.0;
	AxNP[0][1] = 1.0;
	AxNP[0][2] = 0.0;
	AxNP[0][3] = 0.0;
	AxNP[1][0] = -(v1 * v1) + (((gamma - 1.0) * norma_U23 * 0.5) / (rho * rho));
	AxNP[1][1] = (3.0 - gamma) * v1;
	AxNP[1][2] = -(gamma - 1.0) * v2;
	AxNP[1][3] = gamma - 1.0;
	AxNP[2][0] = - v1*v2;
	AxNP[2][1] = v2;
	AxNP[2][2] = v1;
	AxNP[2][3] = 0.0;
	AxNP[3][0] = (-gamma * v1 * E) + ((gamma - 1.0) * v1 * (norma_U23 / (rho * rho)));
	AxNP[3][1] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v1 * v1)));
	AxNP[3][2] = -((gamma - 1.0) * v1 * v2);
	AxNP[3][3] = gamma * v1;

	// *** Ay coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	AyNP[0][0] = 0.0;
	AyNP[0][1] = 0.0;
	AyNP[0][2] = 1.0;
	AyNP[0][3] = 0.0;
	AyNP[1][0] = AxNP[2][0];
	AyNP[1][1] = v2;
	AyNP[1][2] = v1;
	AyNP[1][3] = 0.0;
	AyNP[2][0] = -(v2 * v2) + (((gamma - 1.0) * norma_U23 * 0.5) / (rho * rho)); 
	AyNP[2][1] = -(gamma - 1.0) * v1;
	AyNP[2][2] = (3.0 - gamma) * v2;
	AyNP[2][3] = AxNP[1][3];
	AyNP[3][0] = (- gamma * v2 * E) + (((gamma - 1.0) * v2 * norma_U23) / (rho * rho));
	AyNP[3][1] = AxNP[3][2];
	AyNP[3][2] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v2 * v2)));
	AyNP[3][3] = gamma * v2;
	
	//Produto PAxNP = Ax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for(i=0; i<4; i++)
		for(j=0; j<4; j++){
			Ax[i][j] = 0.0;
			for(k=0; k<4; k++)
				Ax[i][j] += P[i][k] * AxNP[k][j];
		}	
	//Produto PAy = Ay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	for(i=0; i<4; i++)
		for(j=0; j<4; j++){
			Ay[i][j] = 0.0;
			for(k=0; k<4; k++)
				Ay[i][j] += P[i][k] * AyNP[k][j];
		}	
}

void NPreconditioned_Ax_Ay_calculations(double cv, double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4], double Dpmax)
{
	double rho = Ub[0];
	double v1 = Ub[1] / Ub[0];
	double v2 = Ub[2] / Ub[0];
	double E  = Ub[3] / Ub[0];
	double norma_U23 = (Ub[1] * Ub[1]) + (Ub[2] * Ub[2]);

	// *** Ax coefficients
	Ax[0][0] = 0.0;
	Ax[0][1] = 1.0;
	Ax[0][2] = 0.0;
	Ax[0][3] = 0.0;
	Ax[1][0] = -(v1 * v1) + (((gamma - 1) * norma_U23 * 0.5) / (rho * rho));
	Ax[1][1] = (3.0 - gamma) * v1;
	Ax[1][2] = -(gamma - 1.0) * v2;
	Ax[1][3] = gamma - 1.0;
	Ax[2][0] = - v1*v2;
	Ax[2][1] = v2;
	Ax[2][2] = v1;
	Ax[2][3] = 0.0;
	Ax[3][0] = (-gamma * v1 * E) + ((gamma - 1.0) * v1 * (norma_U23 / (rho * rho)));
	Ax[3][1] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v1 * v1)));
	Ax[3][2] = -((gamma - 1.0) * v1 * v2);
	Ax[3][3] = gamma * v1;

	// *** Ay coefficients

	Ay[0][0] = 0.0;
	Ay[0][1] = 0.0;
	Ay[0][2] = 1.0;
	Ay[0][0] = 0.0;
	Ay[1][0] = Ax[2][0];
	Ay[1][1] = v2;
	Ay[1][2] = v1;
	Ay[1][3] = 0.0;
	Ay[2][0] = -(v2 * v2) + (((gamma - 1.0) * norma_U23 * 0.5) / (rho * rho)); 
	Ay[2][1] = -(gamma - 1.0) * v1;
	Ay[2][2] = (3.0 - gamma) * v2;
	Ay[2][3] = Ax[1][3];
	Ay[3][0] = (- gamma * v2 * E) + (((gamma - 1.0) * v2 * norma_U23) / (rho * rho));
	Ay[3][1] = Ax[3][2];
	Ay[3][2] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v2 * v2)));
	Ay[3][3] = gamma * v2;
}	

/*
void dimensionless_Ax_Ay_calculations(double cv, double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4])
{
	// Jacobian matrices Ax and Ay : primitive variables
	double v1 = Ub[1] / Ub[0];
	double v2 = Ub[2] / Ub[0];
	double E  = Ub[3] / Ub[0];
	double norma_v = v1*v1 + v2*v2; //||v||^2

	//Dimensionless Ax and Ay/
	double epsilon, epsilon2;
	epsilon = Mach; //reference Mach number
	epsilon2 = epsilon * epsilon;

	// Ax coefficients
	Ax[0][0] = 0.0;
	Ax[0][1] = 1.0;
	Ax[0][2] = 0.0;
	Ax[0][3] = 0.0;
	Ax[1][0] = -(v1 * v1) + ((gamma - 1) * norma_v * 0.5) ;
	Ax[1][1] = (3.0 - gamma) * v1;
	Ax[1][2] = -(gamma - 1.0) * v2;
	Ax[1][3] = (gamma - 1.0)/epsilon2;
	Ax[2][0] = - v1*v2;
	Ax[2][1] = v2;
	Ax[2][2] = v1;
	Ax[2][3] = 0.0;
	Ax[3][0] = - v1 * (gamma * E - ((gamma - 1.0) * epsilon2 * norma_v));
	Ax[3][1] = gamma * E - (gamma - 1.0) * epsilon2 * ( 0.5 * norma_v + v1 * v1 );
	Ax[3][2] = -(gamma - 1.0) * epsilon2 * v1 * v2;
	Ax[3][3] = gamma * v1;

	// Ay coefficients
	Ay[0][0] = 0.0;
	Ay[0][1] = 0.0;
	Ay[0][2] = 1.0;
	Ay[0][0] = 0.0;
	Ay[1][0] = Ax[2][0];
	Ay[1][1] = v2;
	Ay[1][2] = v1;
	Ay[1][3] = 0.0;
	Ay[2][0] = -(v2 * v2) + (((gamma - 1.0) * norma_v * 0.5)); 
	Ay[2][1] = -(gamma - 1.0) * v1;
	Ay[2][2] = (3.0 - gamma) * v2;
	Ay[2][3] = Ax[1][3];
	Ay[3][0] = - v2 * (gamma * E - (gamma - 1.0) * epsilon2 * norma_v);
	Ay[3][1] = Ax[3][2];
	Ay[3][2] = gamma * E - (gamma - 1.0) * epsilon2 * (0.5 * norma_v + v2 * v2);
	Ay[3][3] = gamma * v2;

}
*/
