function [Diag,invDiag,PA,PF] = monta_Diag(A,F,Elements,Nodes,lm,nodes,neq,nel)

  Diag = zeros(nodes,9);
  invDiag = zeros(nodes,9);
  PA = zeros(nel,81);
  PF = zeros(neq+1,1);
  F(neq+1) = 0;
  PF = F;
  DD = [1 2 3 10 11 12 19 20 21; 
        31	32	33 40	41	42 49	50	51; 
	      61	62	63 70	71	72 79	80	81];
 
  Aux = Nodes;
  
  for e=1:nel
    I1 = Elements(e,1);
    I2 = Elements(e,2);
    I3 = Elements(e,3);
    for i=1:3
      for j=1:3
        Diag(I1,3*(i-1)+j) += A(e,DD(1,3*(i-1)+j));
        Diag(I2,3*(i-1)+j) += A(e,DD(2,3*(i-1)+j));
        Diag(I3,3*(i-1)+j) += A(e,DD(3,3*(i-1)+j));
      end
    end
  end	
  
  for i=1:nodes
      M = reshape(Diag(i,:),3,3)';
      invM = inv(M)';
      invDiag(i,:) = reshape(invM,1,9);
  end
  
  for e=1:nel
      I1 = Elements(e,1);
      I2 = Elements(e,2);
      I3 = Elements(e,3);
      
      Inv11 = (reshape(invDiag(I1,:),3,3)');
      Inv12 = (reshape(invDiag(I1,:),3,3)');
      Inv13 = (reshape(invDiag(I1,:),3,3)');
      Inv21 = (reshape(invDiag(I2,:),3,3)');
      Inv22 = (reshape(invDiag(I2,:),3,3)');
      Inv23 = (reshape(invDiag(I2,:),3,3)');
      Inv31 = (reshape(invDiag(I3,:),3,3)');
      Inv32 = (reshape(invDiag(I3,:),3,3)');
      Inv33 = (reshape(invDiag(I3,:),3,3)');

      Me11 = (reshape(A(e,:),9,9)'(1:3,1:3));
      Me12 = (reshape(A(e,:),9,9)'(1:3,4:6));
      Me13 = (reshape(A(e,:),9,9)'(1:3,7:9));
      Me21 = (reshape(A(e,:),9,9)'(4:6,1:3));
      Me22 = (reshape(A(e,:),9,9)'(4:6,4:6));
      Me23 = (reshape(A(e,:),9,9)'(4:6,7:9));
      Me31 = (reshape(A(e,:),9,9)'(7:9,1:3));
      Me32 = (reshape(A(e,:),9,9)'(7:9,4:6));
      Me33 = (reshape(A(e,:),9,9)'(7:9,7:9));

      if (Aux(I1,1)>neq)
	      Me11(1,:) = [1 0 0];
	      Me11(:,1) = [1 0 0]';
	      Me12(1,:) = [0 0 0];
	      Me13(1,:) = [0 0 0];	
      end
     if (Aux(I1,2)>neq)
      	Me11(2,:) = [0 1 0];
	      Me11(:,2) = [0 1 0]';
	      Me12(2,:) = [0 0 0];
	      Me13(2,:) = [0 0 0];	
     end
     if (Aux(I1,3)>neq)
	      Me11(3,:) = [0 0 1];
	      Me11(:,3) = [0 0 1]';
	      Me12(3,:) = [0 0 0];
	      Me13(3,:) = [0 0 0];	
      end
     if (Aux(I2,1)>neq)
	      Me21(1,:) = [0 0 0];
      	Me22(1,:) = [1 0 0];
	      Me22(:,1) = [1 0 0]';
	      Me23(1,:) = [0 0 0];	
     end
     if (Aux(I2,2)>neq)
      	Me21(2,:) = [0 0 0];
	      Me22(2,:) = [0 1 0];
	      Me22(:,2) = [0 1 0]';
	      Me23(2,:) = [0 0 0];	
      end
     if (Aux(I2,3)>neq)
	      Me21(3,:) = [0 0 0];
	      Me22(3,:) = [0 0 1];
	      Me22(:,3) = [0 0 1]';
	      Me23(3,:) = [0 0 0];	
     end
       	
     if (Aux(I3,1)>neq)
	      Me31(1,:) = [0 0 0];
	      Me32(1,:) = [0 0 0];	
	      Me33(1,:) = [1 0 0];
	      Me33(:,1) = [1 0 0]';
     end
     
     if (Aux(I3,2)>neq)
	      Me31(2,:) = [0 0 0];
	      Me32(2,:) = [0 0 0];	
	      Me33(2,:) = [0 1 0];
	      Me33(:,2) = [0 1 0]';
      end
     if (Aux(I3,3)>neq)
	      Me31(3,:) = [0 0 0];
	      Me23(3,:) = [0 0 0];	
	      Me33(3,:) = [0 0 1];
	      Me33(:,3) = [0 0 1]';
     end
     
      PAe = [(Inv11*Me11) (Inv12*Me12) (Inv13*Me13); (Inv21*Me21) (Inv22*Me22) (Inv23*Me23); (Inv31*Me31) (Inv32*Me32) (Inv33*Me33)]';
      PA(e,:) = reshape(PAe,1,81);
      
   end

  for i=1:nodes
      
      a = F(Aux(i,1));
      b = F(Aux(i,2));
      c = F(Aux(i,3));
      
      F(Aux(i,1)) = invDiag(i,1:3)* ([a b c]');	
      F(Aux(i,2)) = invDiag(i,4:6)* ([a b c]');	
      F(Aux(i,3)) = invDiag(i,7:9)* ([a b c]');	
            
      F(neq+1) = 0;
      
  end
   
  

  PF = F(1:neq);

