function [Xi, B, Psi] = process_eigenfunction(V, index, Kx, Kz, Fr,...
    alpha, beta, Lx, Lz)

%% SELECT DESIRED VECTOR

qhat = V(:, index);

%% DECOMPOSE EIGENFUNCTION

[xihat, bhat] = decompose_eigenfunction(qhat); 

%% GET Nx, Nz 

Nx = size(Kx,2);
Nz = size(Kx,1);

%% CALCULATE psihat

psihat = inverse_laplacian(xihat, Kx, Kz, Fr, alpha, beta, Lx, Lz); 

%% RESHAPE EIGENVECTORS

Xihat2  = reshape(xihat, Nz, Nx); 
Bhat2   = reshape(bhat, Nz, Nx);  
Psihat2 = reshape(psihat, Nz, Nx);

%% IFFT

Xi  = real(ifft2(Xihat2));
Psi = real(ifft2(Psihat2));
B   = real(ifft2(Bhat2));

end