# include "naca0012.h"

int NACA0012_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){

	int I, nnodes;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++){
		
		for(I = 0; I < nnodes; I++)
	{
		if (Node[I].id[0] >= 0){ u[Node[I].id[0]] = 1.0; }
		if (Node[I].id[1] >= 0){ u[Node[I].id[1]] = 1.0; }
		if (Node[I].id[2] >= 0){ u[Node[I].id[2]] = 0.0; }
		if (Node[I].id[3] >= 0){ u[Node[I].id[3]] = 10.5; } // c_v = 10 para definir mach = 0.5
	}

	}

	return 0;
}


