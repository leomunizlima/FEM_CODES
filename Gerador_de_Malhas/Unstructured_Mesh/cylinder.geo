//parametro de malha
//lc = 0.7e-1;
lc = 5e-2; 
Point(1) = {-3, 0, 0, lc};
Point(2) = {-0.5, 0, 0, lc};
Point(3) = {0, 0, 0, lc};
Point(4) = {0, 0.5, 0, lc};
Point(5) = {0.5, 0, 0, lc};

Point(6) = {3, 0, 0, lc};
Point(7) = {3, 3, 0, lc};
Point(8) = {-3, 3, 0, lc};

Line(1) = {1, 2};
Circle(2) = {2, 3, 4};
Circle(3) = {4, 3, 5};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 1};

Line Loop(8) = {1, 2, 3, 4, 5, 6, 7}; 
Plane Surface(9) ={8}; 

//Plane Surface(10) ={9}; 
//Plane Surface(9) = {7}; 

//Physical Point(1) = {1,2};

//Dirichlet = 1;				//marca dos nós de fronteira
//Neumann = 2;					//marca dos nós de fronteira
//Physical Line(Dirichlet) = {2,5,6};		//nós de fronteira
//Physical Line(Neumann) = {1,3,4}; 		//nós de fronteira
//Physical Point(Dirichlet) = {1,4,5,6};	//nós de fronteira
//Physical Point(Neumann) = {2,3};		//nós de fronteira
//Physical Surface("My fancy surface label") = {9}
Physical Surface("My surface") = {9} ;

//Color Grey50{ Surface{ 8 }; }
//Color Purple{ Surface{ 9 }; }
//Color Red{ Line{ 1:4 }; }
//Color Yellow{ Line{ 5:6 }; }

