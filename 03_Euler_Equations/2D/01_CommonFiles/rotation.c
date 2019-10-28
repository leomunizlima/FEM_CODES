#include "EulerEquations.h"

void rotation(int tag, double theta, double M[12][12], double K[12][12])
{
	int I, J, k, cont1, cont2;
	double sum1, sum2, R[12][12], RT[12][12];
	
	cont1 = 4*tag-4;
	cont2 = 4*tag;
	
	R[cont1+0][cont1+0] = 1.0;
	R[cont1+0][cont1+1] = 0.0;
	R[cont1+0][cont1+2] = 0.0;
	R[cont1+0][cont1+3] = 0.0;

	R[cont1+1][cont1+0] = 0.0;
	R[cont1+1][cont1+1] = cos(theta);
	R[cont1+1][cont1+2] = -sin(theta);
	R[cont1+1][cont1+3] = 0.0;

	R[cont1+2][cont1+0] = 0.0;
	R[cont1+2][cont1+1] = sin(theta);
	R[cont1+2][cont1+2] = cos(theta);
	R[cont1+2][cont1+3] = 0.0;

	R[cont1+3][cont1+0] = 0.0;
	R[cont1+3][cont1+1] = 0.0;
	R[cont1+3][cont1+2] = 0.0;
	R[cont1+3][cont1+3] = 1.0;
	
	//TRANSPOSTA
	RT[cont1+0][cont1+0] = 1.0;
	RT[cont1+0][cont1+1] = 0.0;
	RT[cont1+0][cont1+2] = 0.0;
	RT[cont1+0][cont1+3] = 0.0;

	RT[cont1+1][cont1+0] = 0.0;
	RT[cont1+1][cont1+1] = cos(theta);
	RT[cont1+1][cont1+2] = sin(theta);
	RT[cont1+1][cont1+3] = 0.0;

	RT[cont1+2][cont1+0] = 0.0;
	RT[cont1+2][cont1+1] = -sin(theta);
	RT[cont1+2][cont1+2] = cos(theta);
	RT[cont1+2][cont1+3] = 0.0;

	RT[cont1+3][cont1+0] = 0.0;
	RT[cont1+3][cont1+1] = 0.0;
	RT[cont1+3][cont1+2] = 0.0;
	RT[cont1+3][cont1+3] = 1.0;

	//A*R
	for(I=0; I<12 ;I++)
		for(J=cont1; J<cont2; J++){
			sum1 = 0.0;
			sum2 = 0.0;
			for(k=cont1; k<cont2; k++){
				sum1 +=M[I][k]*R[k][J]; 	
				sum2 +=K[I][k]*R[k][J];
			}
			M[I][J] = sum1;
			K[I][J] = sum2;
		}
	//RT*A
	for(I=cont1; I<cont2; I++)
		for(J=0; J<12 ;J++){
			sum1 = 0.0;
			sum2 = 0.0;
			for(k=cont1; k<cont2; k++){
				sum1 +=RT[I][k]*M[k][J]; 	
				sum2 +=RT[I][k]*K[k][J];
			}
			M[I][J] = sum1;
			K[I][J] = sum2;
		}
	/*
	for(I=0; I<12 ;I++)
		for(J=cont1; J<cont2; J++){
			sum1 = 0.0;
			sum2 = 0.0;
			for(k=0; k<12; k++){
				sum1 +=RT[I][k]*M[k][J]; 	
				sum2 +=RT[I][k]*K[k][J];
			}
			M[I][J] = sum1;
			K[I][J] = sum2;
		}
	*/
	//for(I=0; I<12 ;I++){
	//	sum1 = 0.0;
	//	for(k=0; k<12; k++)
	//		sum1 += RT[I][k]*RES[k];			
	//	RES[I] = sum1;		
	//}
	return;	

}




