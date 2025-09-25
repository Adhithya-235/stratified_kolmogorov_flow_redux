function [XiT, BT, PsiT] = process_eigenfunction_2(V, index, Kx, Kz, Fr,...
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

Xi  = ifft2(Xihat2);
Psi = ifft2(Psihat2);
B   = ifft2(Bhat2);

%% CONCATENATE

XiC  = cat(2, Xi, Xi);
PsiC = cat(2, Psi, Psi);
BC   = cat(2, B, B);

%% COMPUTE x AND z, FLOQUET MULTIPLIER

dx    = Lx/Nx;
dz    = Lz/Nz;
x     = 0:dx:(2*Lx)-dx;
z     = 0:dz:Lz-dz;
[X,Z] = meshgrid(x,z);
k     = 2*pi/Lx;
l     = 2*pi/Lz;
floqe = exp(1i*(alpha*k*X + beta*l*Z));

%% CONSTRUCT FLOQUET MODIFIED PERTURBATION

XiL  = floqe.*XiC;
PsiL = floqe.*PsiC;
BL   = floqe.*BC;

%% INTERPOLATE TO NATIVE SIZE

XiT  = interpft(interpft(real(XiL), Nz, 1), Nx, 2); 
PsiT = interpft(interpft(real(PsiL), Nz, 1), Nx, 2); 
BT   = interpft(interpft(real(BL), Nz, 1), Nx, 2); 

end