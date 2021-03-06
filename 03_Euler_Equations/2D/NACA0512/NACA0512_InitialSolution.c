# include "naca0512.h"

int NACA0512_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u)
{

	int I, nnodes;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].id[0] >= 0){ 
			u[Node[I].id[0]] = 1.0; 
		}
		if (Node[I].id[1] >= 0){ 
			//if (Node[I].v1Type<0){
			//	theta = NACA0512_theta(Node[I].x, Node[I].y);
			//	u[Node[I].id[1]] = cos(theta); //ajustando a condição inicial para a malha do naca			
			//}
			//else	
				u[Node[I].id[1]] = cos(0.0872664625997165);
		}
		if (Node[I].id[2] >= 0){ 
			u[Node[I].id[2]] = sin(0.0872664625997165);
		}
		if (Node[I].id[3] >= 0){ 
			u[Node[I].id[3]] = NACA0512_epresc(Parameters,0,0); // we used x=0 and y=0 because x and y are not necessary
		}	

/*		if (Node[I].id[3] >= 0){ 
			u[Node[I].id[3]] = 1785714.7857143;  // rho * e = 0.9464 -> Mach=2 (Cv = 716.5)
		}	//rho * e = 2.28571->mach = 1.0 (Cv = 1.78571) %%% rho * e = 2.78795->mach = 0.8 (Cv = 2.28795) %%% rho * e = 7.642857->mach = 0.5 (Cv = 7.142857) %%% rho * e = 179.07143->mach = 0.1 (Cv=178.57143) 
			//%%% rho * e = 20.34127->mach = 0.3 (Cv=19.84127) %%% rho * e = 17857.64286->mach = 0.01 (Cv=17857.14286)
	}
*/
	}	
	return 0;
}


