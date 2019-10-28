#include "EulerEquations.h"

void rotationNMV(int tag, double theta, double Me[12][12], double Mhh[12][12], double Khh[12][12], double KhB[12][4], double N1[12][4] , double MBh[4][12], double KBh[4][12]){
	int I, J, k, cont1, cont2;
	double sum1, sum2, sum3, R[12][12], RT[12][12];
	
	cont1 = 4*(tag-1);
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

	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	//Me*R, Mhh*R e Khh*R
	for(I=0; I<12 ;I++)
		for(J=cont1; J<cont2; J++){
			sum1 = 0.0;	
			sum2 = 0.0;	
			sum3 = 0.0;
			for(k=cont1; k<cont2; k++){
				sum1 +=Me[I][k]*R[k][J];	
				sum2 +=Mhh[I][k]*R[k][J];		
				sum3 +=Khh[I][k]*R[k][J];
			}
			Me[I][J] = sum1;			
			Mhh[I][J] = sum2;
			Khh[I][J] = sum3;
		}
	//RT*Me, RT*Mhh	e RT*Khh
	for(I=cont1; I<cont2; I++)
		for(J=0; J<12 ;J++){
			sum1 = 0.0;
			sum2 = 0.0;
			sum3 = 0.0;
			for(k=cont1; k<cont2; k++){
				sum1 +=RT[I][k]*Me[k][J]; 	
				sum2 +=RT[I][k]*Mhh[k][J];
				sum3 +=RT[I][k]*Khh[k][J];
			}
			Me[I][J] = sum1;
			Mhh[I][J] = sum2;
			Khh[I][J] = sum3;
		}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	//MBh*R, KBh*R
	for(I=0; I<4 ;I++)
		for(J=cont1; J<cont2; J++){						
			sum2 = 0.0;
			sum3 = 0.0;
			for(k=cont1; k<cont2; k++){									
				sum2 +=MBh[I][k]*R[k][J];
				sum3 +=KBh[I][k]*R[k][J];
			}						
			MBh[I][J] = sum2;
			KBh[I][J] = sum3;
		}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	//RT*KhB e RT*N1
	for(I=cont1; I<cont2 ;I++)
		for(J=0; J<4; J++){						
			sum1 = 0.0;			
			sum2 = 0.0;
			for(k=cont1; k<cont2; k++){									
				sum1 +=RT[I][k]*KhB[k][J];
				sum2 +=RT[I][k]*N1[k][J];
			}						
			KhB[I][J] = sum1;			
			N1[I][J] = sum2;
		}
	
	return;	

}




