function [M] = monta_Matriz_Global(A,lm, neq, nel)
	M = sparse(neq,neq);
	for e=1:nel
		for i=1:9
			for j=1:9
				if (lm(e,i)<=neq && lm(e,j)<=neq)
					M(lm(e,i),lm(e,j)) += A(e,9*(i-1)+j);
				end
			end
		end
	end
	
end



