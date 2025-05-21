function [xihat, bhat] = decompose_eigenfunction(qhat)

MN    = length(qhat)/2;
xihat = qhat(1:MN);
bhat  = qhat((MN+1):(2*MN));

end