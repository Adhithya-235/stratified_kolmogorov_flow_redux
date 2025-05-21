function [Xi, B, Psi] = process_eigenfunction(V, index, Kx, Kz, Fr,...
    alpha, beta, Lx, Lz)

%% SELECT DESIRED VECTOR

q = V(:, index);

%% DECOMPOSE EIGENFUNCTION

[xi, b] = decompose_eigenfunction(q); 

%% CALCULATE psihat

psi = inverse_laplacian(xi, Kx, Kz, Fr, alpha, beta, Lx, Lz);

%% GET Nx, Nz 

Nx = size(Kx,2);
Nz = size(Kx,1);

%% RESHAPE EIGENVECTORS

Xi  = reshape(xi, Nz, Nx); 
B   = reshape(b, Nz, Nx);  
Psi = reshape(psi, Nz, Nx);

end