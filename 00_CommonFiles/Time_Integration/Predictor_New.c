#include "time_integration.h"

//int Special_Paraview_Output(int, ParametersType *, FemStructsType *, FemFunctionsType *, double *, double *, double *, double *);

int Predictor_New(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int I, i;
	int nel, neq, passo;
	double t, dt, alpha, norm_a, norm_Da, tol_correction;
	double *a, *aB, *auxVec, *Da, *DaB, *u, *u_old, **R2, *R2aux, *invN2, **M2, *M2aux, *F; //Parametros do Preditor
	double *uB, *delta_old, *deltaNMV_old;
	double **R1, *R1aux;
	AuxBuildStructuresType *AuxBuild;
	
	nel = Parameters->nel;
	neq = Parameters->neq;
	
	a = (double*) mycalloc("a of 'Preditor_New'", neq + 1, sizeof(double));
	aB = (double*) mycalloc("aB of 'Preditor_New'",NDOF*nel, sizeof(double));
	Da = (double*) mycalloc("Da of 'Preditor_New'", neq + 1, sizeof(double));
	DaB = (double*) mycalloc("DaB of 'Preditor_New'", NDOF*nel, sizeof(double));
	auxVec = (double*) mycalloc("auxVec of 'Preditor_New'", neq + 1, sizeof(double));
	u_old = (double*) mycalloc("u_old of 'Preditor_New'", neq + 1, sizeof(double));
	R2 = (double**) mycalloc("R2 of 'Preditor_New'", nel, sizeof(double));
	R2aux = (double *) mycalloc("R2aux of 'Preditor_New'", NDOF*nel, sizeof(double));
	for (I = 0; I < nel;I++)
		R2[I] = &R2aux[NDOF*I]; 
	R1 = (double**) mycalloc("R1 of 'Preditor_New'", nel, sizeof(double));
	R1aux = (double *) mycalloc("R1aux of 'Preditor_New'", 12*nel, sizeof(double));
	for (I = 0; I < nel;I++)
		R1[I] = &R1aux[12*I];
	invN2 = (double*) mycalloc("invN2 of 'Preditor_New'", nel, sizeof(double));
	M2 = (double **) mycalloc("M2 of 'Preditor_New'", nel, sizeof(double*));
	M2aux = (double *) mycalloc("M2aux of 'Preditor_New'", NNOEL*NDOF*NDOF*nel, sizeof(double));
	for (I = 0; I < nel;I++)
		M2[I] = &M2aux[NNOEL*NDOF*NDOF*I]; 
	uB = (double*) mycalloc("uB of 'Preditor_New'", nel*NDOF, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Preditor_New'", nel, sizeof(double));
	deltaNMV_old = (double*) mycalloc("deltaNMV_old of 'Preditor_New'", nel, sizeof(double));

	u = FemStructs->u;
	F = FemStructs->F;
	dt = Parameters->DeltaT;
	Parameters->DeltaT_Build = dt;
	alpha = Parameters->Alpha;
	Parameters->Alpha_Build = alpha;
	tol_correction = Parameters->NonLinearTolerance;
	AuxBuild = (AuxBuildStructuresType*) mycalloc("AuxBuild of 'Predictor_New'",1,sizeof(AuxBuildStructuresType));
	AuxBuild->tolerance = Parameters->StabilizationTolerance;
	AuxBuild->M2 = M2;
	AuxBuild->R1 = R1;
	AuxBuild->R2 = R2;
	AuxBuild->invN2 = invN2;
	FemStructs->AuxBuild = AuxBuild;
	FemStructs->delta_old = delta_old;
	FemStructs->deltaNMV_old = deltaNMV_old;
	FemStructs->du = a;
	FemStructs->uB = uB;
	FemStructs->duB = aB;

	FemFunctions->InitialSolution(Parameters, FemStructs->Node, u);
	#ifdef EulerEquations2D
		int nnodes=Parameters->nnodes;
		double *U, *U1, *U2, *rho, *e, *diff;
		double x,y;

		double res_med = 0.0;
		double res_min = 1e20;
		double res_max = 0;
		double res;

		int n_passos = 30;
	//	int count = 0;

		double v_min = 0.8;
		double v_max = 1.2;


		U = (double*) mycalloc("U1 of 'Predictor_New'",nnodes,sizeof(double));
		U1 = (double*) mycalloc("U1 of 'Predictor_New'",nnodes,sizeof(double));
		U2 = (double*) mycalloc("U1 of 'Predictor_New'",nnodes,sizeof(double));
		rho = (double*) mycalloc("rho of 'Predictor_New'",nnodes,sizeof(double));
		e = (double*) mycalloc("e of 'Predictor_New'",nnodes,sizeof(double));
		diff = (double*) mycalloc("diff of 'Predictor_New'",nnodes,sizeof(double));
	
		 
		for (I=0;I<nnodes;I++){
		
				x = FemStructs->Node[I].x;
				y = FemStructs->Node[I].y;
			
				if (FemStructs->Node[I].id[0]==-1)
					rho[I] = FemFunctions->rhopresc(x,y);
			
				if (FemStructs->Node[I].id[1]==-1)
					U1[I] = FemFunctions->v1presc(x,y);
			
				if (FemStructs->Node[I].id[2]==-1)
					U2[I] = FemFunctions->v2presc(x,y);
			
				if (FemStructs->Node[I].id[3]==-1)
					e[I]= FemFunctions->epresc(Parameters,x,y);
			}
		#endif

		t = 0.0;
		passo = 0;
		int tag = 1;
		
		do{
			passo++;
			t += dt;
			
			#ifdef debug
				printf("\n\n Passo: %d\n", passo); 
			#endif
			//PREDICAO
			i = 0;
			
			dmemcpy(neq, u, u_old); // copy u to u_old		
			daxpy(neq, (1.0-alpha)*dt , a, u); // u  = u + (1-alpha)*dt*a
			memsetzero(neq,a); // set 0 in a
		
			daxpy(nel*NDOF, (1.0-alpha)*dt , aB, uB); // uB  = uB + (1-alpha)*dt*aB
			memsetzero(nel*NDOF,aB); // set 0 in aB
		
			// MULTICORRECAO
			do{
				i++;

				FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
			
				FemFunctions->scaling(Parameters, MatrixData, FemStructs);

				FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, F);

				FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, Da);
				
				FemFunctions->unscaling(Parameters, MatrixData, FemStructs, Da);

				calculate_DaB(Parameters, FemStructs, FemFunctions, Da, DaB);

				daxpy(neq, 1, Da, a);
				daxpy(neq, alpha*dt, Da, u);

				daxpy(nel*NDOF, 1, DaB, aB);
				daxpy(nel*NDOF, alpha*dt, DaB, uB);

				norm_a = sqrt(ddot(neq, a, a));
				norm_Da = sqrt(ddot(neq, Da, Da));
				#ifdef debug
					double normF;
					normF = sqrt(ddot(neq, F, F));
					printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n", 
						tol_correction*norm_a, normF, norm_a, norm_Da, t, i);
				#endif
			}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection
			
			
			#ifdef debug
				printf("\n\n"); 
			#endif

						
			#ifdef EulerEquations2D

			double norm2_rho, norm2_U, norm2_e, norm2_diff;  
			for (I=0;I<nnodes;I++){
				if (FemStructs->Node[I].id[0]!=-1){
					rho[I] = u[FemStructs->Node[I].id[0]];
					diff[I] = fabs(rho[I] - u_old[FemStructs->Node[I].id[0]]);
				}	
				else
					diff[I] = 0;
					
			}
				
			norm2_rho = sqrt(ddot(nnodes,rho,rho));
			norm2_diff = sqrt(ddot(nnodes,diff,diff));
			norm2_rho = norm2_diff/norm2_rho;
				
			double diffU1, diffU2;
			for (I=0;I<nnodes;I++){
				if (FemStructs->Node[I].id[1]!=-1){
					U1[I] = u[FemStructs->Node[I].id[1]];
					diffU1 = fabs(U1[I]-u_old[FemStructs->Node[I].id[1]]);
				}	
				else
					diffU1 = 0;

				if (FemStructs->Node[I].id[2]!=-1){
					U2[I] = u[FemStructs->Node[I].id[2]];
					diffU2 = fabs(U2[I]-u_old[FemStructs->Node[I].id[2]]);
				}	
				else
					diffU2 = 0;

				U[I] = sqrt(U1[I]*U1[I]+U2[I]*U2[I]);
				diff[I] = sqrt(diffU1*diffU1 + diffU2*diffU2);

			}
			
			norm2_U = sqrt(ddot(nnodes,U,U));
			norm2_diff = sqrt(ddot(nnodes,diff,diff));
			norm2_U = norm2_diff/norm2_U;

			for (I=0;I<nnodes;I++){
				if (FemStructs->Node[I].id[3]!=-1){
					e[I] = u[FemStructs->Node[I].id[3]];
					diff[I] = fabs(e[I]-u_old[FemStructs->Node[I].id[3]]);
				}	
				else
					diff[I] = 0;
				
			}
			
			norm2_e = sqrt(ddot(nnodes,e,e));
			norm2_diff = sqrt(ddot(nnodes,diff,diff));
			norm2_e = norm2_diff/norm2_e;

			fprintf(stderr,"%.3lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n",t,norm2_U,norm2_rho,norm2_e, norm2_U + norm2_rho + norm2_e);
		
			if (norm2_e < 1e-7) break; 

			//Residue freezing


			res = norm2_rho;

			res_med += res; 

			if (res < res_min)
				res_min = res;

			if (res > res_max)
				res_max = res;	

			if (passo % n_passos == 0){
				res_med = res_med / n_passos;


				if (res_min >= v_min * res_med && res_max <= v_max * res_med){
					FemFunctions->ShockCapture = FemFunctions->ShockCapture_NOT;
		//			FemFunctions->ShockCapture_Micro = FemFunctions->ShockCapture_NOT;
				//	printf("\n\n\n PASSO=%d res_med=%lf  res_min=%lf    res_max=%lf\n\n\n",passo,res_med,res_min,res_max);
				//	count++;
				/*	if (count>30){
						setDelta(Parameters, FemFunctions, FemOtherFunctions);
						count = 0;
					}	*/

				}
		//		else
		//			setDelta(Parameters, FemFunctions, FemOtherFunctions);
				res_med = 0;
				res_min = 1e20;
				res_max = 0;

				Special_Paraview_Output(passo, Parameters, FemStructs, FemFunctions, rho, U1, U2, e);
			}	
			

		#endif	



	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time

	free(a);
	free(Da);
	free(u_old);
	free(auxVec);
	free(uB);
	free(DaB);
	free(delta_old);
	free(M2aux);
	free(M2);
	free(R2aux);
	free(R2);
	free(R1aux);
	free(R1);
	free(invN2);
	free(deltaNMV_old);
	free(AuxBuild);


	#ifdef EulerEquations2D
	free(rho);
	free(U);
	free(U1);
	free(U2);
	free(e);
	free(diff);
	#endif

	return 0;

}


/*
int Special_Paraview_Output(int passo, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *rho, double *U1, double *U2, double *e)
{
	int I, nnodes, nel;
	char FileName[2000];
	FILE *OutFile;
	double X, Y, aux;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	nnodes = Parameters->nnodes;
	nel = Parameters->nel;
	double theta,gamma;	//my alteration
	double *pres;

	pres = (double*) mycalloc("pres of 'Paraview_Output'", nnodes, sizeof(double));

	for (I = 0; I < nnodes; I++){
		X = Node[I].x;
		Y = Node[I].y;
		gamma = FemFunctions->gamma(X, Y);

		if (Node[I].id[1] >= 0 && Node[I].v1Type < 0){				//my alteration
			theta = FemFunctions->BC_theta(X,Y);
			U1[I] = cos( theta ) * aux;	//projection Us -> v1 = cos(theta) * Us	
		}
	   	
	   	if (Node[I].id[2] >= 0 && Node[I].v1Type < 0){
			theta = FemFunctions->BC_theta(X,Y);
			U2[I] =  sin( theta ) * aux;		//projection Us -> v2 = sin(theta) * Us	
		}
	
		if (Node[I].v1Type < 0){				//my alteration
			theta = FemFunctions->BC_theta(X,Y);
			aux = U1[I];
			U1[I] = cos( theta ) * aux;	//projection Us -> v1 = cos(theta) * Us	
			U2[I] =  sin( theta ) * aux;		//projection Us -> v2 = sin(theta) * Us	
			
		}	
		aux = e[I] - (U1[I]*U1[I] + U2[I]*U2[I])/(2.0*rho[I]);
		pres[I] = (gamma - 1)*aux;
	}

	sprintf(FileName,"/media/leonardo/Arquivos/Linux/VTK/%s_%s_%s_%05d.vtk", Parameters->Experiments, Parameters->ProblemTitle,Parameters->Preconditioned, passo/30);	
	OutFile = myfopen(FileName,"w");

	fprintf(OutFile,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	fprintf(OutFile,"\t<UnstructuredGrid>\n");
	fprintf(OutFile,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nnodes, nel);
	fprintf(OutFile,"\t\t\t<PointData Scalars=\"scalars\" Vectors = \"Velocity\">\n");	

	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", rho[I]);
	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Pression\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", pres[I]);
	fprintf(OutFile,"\t\t\t\t</DataArray>\n");

	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Energy\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", e[I]);
	fprintf(OutFile,"\t\t\t\t</DataArray>\n");

	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\t%.12lf\t%.12lf\n", U1[I], U2[I], 0.0);
	fprintf(OutFile,"\t\t\t\t</DataArray>\n");

	fprintf(OutFile,"\t\t\t</PointData>\n");
	fprintf(OutFile,"\t\t\t<Points>\n");

	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\t%.12lf\t%.12lf\n", Node[I].x,Node[I].y, 0.0);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Points>\n");
	fprintf(OutFile,"\t\t\t<Cells>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
		fprintf(OutFile,"\t\t\t\t   %d\t%d\t%d\n", Element[I].Vertex[0], Element[I].Vertex[1], Element[I].Vertex[2]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
		fprintf(OutFile,"\t\t\t\t   %d\n",(I+1)*3);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
			fprintf(OutFile,"\t\t\t\t   %d\n",5);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Cells>\n");
	fprintf(OutFile,"\t\t</Piece>\n");
	fprintf(OutFile,"\t</UnstructuredGrid>\n");
	fprintf(OutFile,"</VTKFile>\n");

	fclose(OutFile);
	free(pres);

	return 0;
}
*/
