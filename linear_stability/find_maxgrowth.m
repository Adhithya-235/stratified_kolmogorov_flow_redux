function [sigm, vort, buoy, strm] = find_maxgrowth(L, B)

%   Solve the generalized eigenvalue problem LV = sig*BV, remove infinity 
%   eigenvalues, sort the remaining, and return them together with the 
%   corresponding eigenvectors. 

%% SOLVE GENERALIZED EIGENVALUE PROBLEM

[V, sig] = eig(L, B); 

%% DIAGONALIZE EIGENVALUES, SORT BY DESCENDING REAL PART

dsig        = diag(sig);
[~, sindex] = sort(real(dsig), 'descend');
sig_sorted  = dsig(sindex);
V_sorted    = V(:, sindex);

%% REMOVE INFINITIES

sigm = sig_sorted; sigm(isinf(sig_sorted)) = [];
Vm   = V_sorted;   Vm(:,isinf(sig_sorted)) = [];

%% GET VORTICITY, BUOYANCY, STREAMFUNCION

N    = size(L,1)/3;
vort = Vm(1:N,:);
buoy = Vm((N+1):(2*N),:);
strm = Vm(((2*N)+1):(3*N),:);

end