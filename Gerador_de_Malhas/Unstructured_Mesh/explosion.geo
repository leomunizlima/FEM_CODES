//parametro de malha
lc = 1e-1;
Point(1) = {0, 0, 0, lc};
Point(2) = {2, 0, 0, lc};
Point(3) = {2, 2, 0, lc};
Point(4) = {0, 2, 0, lc};
Point(5) = {1, 1.4, 0, lc};
Point(6) = {1, 0.6, 0, lc};
Point(7) = {1, 1, 0, lc};

Line(1) = {2,1};
Line(2) = {1,4};
Line(3) = {4,3};
Line(4) = {3,2};

Circle(5) = {5,7,6};
Circle(6) = {6,7,5};
Line Loop(7) = {5,6};//buraco
Plane Surface(8) ={7}; 

Line Loop(9) = {1,2,3,4}; //retângulo
Plane Surface(10) ={9,7}; 

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

