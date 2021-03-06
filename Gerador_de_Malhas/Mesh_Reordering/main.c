#include "Mesh_Reordering.h"

int main(int argc, char **argv)
{
	int I;
	char Problem[300];
	int nnodes, nel, NDOF;
	ElementType *Element;

	/*****************************************************************************************************************************************/
	//							Testing initial parameters
	/******************************************************************************************************************************************/
	if (argc != 3){
		printf("Use ./Mesh_reordering <file.dat> <NDOF>\n");
		exit(1);
	}

	/*****************************************************************************/
	//				Reading nodes 
	/*****************************************************************************/		
	sscanf(argv[1],"%[^_]_%d_%d.dat",Problem,&nnodes,&nel);

	NDOF = atoi(argv[2]);
	double *X = mycalloc("X",nnodes,sizeof(double));
	double *Y = mycalloc("Y",nnodes,sizeof(double));
	int *TypeAux = mycalloc("TypeAux",nnodes*NDOF,sizeof(int));
	int **Type = mycalloc("Type",nnodes,sizeof(int*));
	for (I=0; I<nnodes; I++){
		Type[I] = &TypeAux[I*NDOF];
	}	
	Element = mycalloc("Element",nel,sizeof(ElementType));

	/************************************************************************************************************/	
	//					Reading mesh
	/************************************************************************************************************/	
	reading_mesh(argv[1],nnodes,nel,NDOF,X,Y,Type,Element);


	/*************************************************************************************************************/
	//					Checking for missing or duplicated points
	/*************************************************************************************************************/	
	checking_mesh(&nnodes, nel, X, Y, Type, Element);


	/**********************************************************************************************************/
	//					Reordering using triangular RCM					
	/**********************************************************************************************************/
	reordering_rcm(nnodes, nel, NDOF, Problem, X, Y, Type, Element);
		
		
	/************************************************************************************************************/	
	//					Writing reordered mesh file
	/************************************************************************************************************/	
	writing_mesh(nnodes, nel, NDOF, Type, Element, X, Y, Problem);


	/************************************************************************************************************/	
	//						Deallocations
	/************************************************************************************************************/	
	free(X);
	free(Y);
	free(TypeAux);
	free(Type);
	free(Element);

	return 0;
}


