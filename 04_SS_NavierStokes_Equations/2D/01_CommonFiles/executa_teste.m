[Diag,invDiag,PA,PF] = monta_Diag(A,F,Elements,Nodes,lm,nodes,neq,nel);
MA = monta_Matriz_Global(A,lm, neq, nel);
MPA = monta_Matriz_Global(PA,lm, neq, nel);
