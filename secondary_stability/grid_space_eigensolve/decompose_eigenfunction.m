function [xi, b] = decompose_eigenfunction(q)

MN = length(q)/2;
xi = q(1:MN);
b  = q((MN+1):(2*MN));

end