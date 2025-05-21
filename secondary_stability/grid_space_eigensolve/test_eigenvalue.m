clear
close all
clc

%% PHYSICAL PARAMETERS

Reb   = 1;                  % Buoyancy Reynolds
Pr    = 1;                  % Prandtl
Fr    = 0.02;               % Froude
Lx    = 0.06*2;             % x-period, base state
Lz    = 4*pi/3;             % z-period, base state
alpha = 0;                  % x Floquet Modifier
beta  = 0;                  % z Floquet Modifier
m     = 3;

%% NUMERICAL PARAMETERS

Nx = 128;
Nz = 128; 

%% SPATIAL GRID

x = linspace(0, Lx, Nx+1); x = x(1:end-1);
z = linspace(0, Lz, Nz+1); z = z(1:end-1);
[X,Z] = meshgrid(x, z);

%% FOURIER GRID

kxVals = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):-1];
kzVals = (2*pi/Lz)*[0:(Nz/2-1) (-Nz/2):-1];
[Kx, Kz] = meshgrid(kxVals, kzVals);

%% GET BASIC STATE

[Xi, B, Psi] = laminar_bs(Z, m);

%% MATRIX OPERATOR AS A FUNCTION HANDLE

Afun = @(q) boussinesq(q, Xi, B, Psi, Kx, Kz, alpha, beta,...
    Lx, Lz, Fr, Pr, Reb);

%% EIGENVALUE SOLVER OPTIONS

dim        = 2*Nx*Nz;
num        = 1;    % COMPUTE 10 LARGEST EIGENVALUES (BY REAL PART)
type       = 'largestreal';
opts.tol   = 1e-12;
opts.maxit = 1000;
opts.disp  = 1;
opts.p     = 20;
opts.v0    = construct_StartVector(Nx, Nz); 

%% EIGENVALUE SOLVE

[V, D] = eigs(Afun, dim, num, type, opts);
disp('Computed eigenvalues:');
disp(Fr*diag(D))

%% PROCESS AND PLOT DOMINANT EIGENFUNCTION

index = 1;
[Xi, B, Psi] = process_eigenfunction(V, index, Kx, Kz, Fr, alpha, beta,...
    Lx, Lz);
% xi2 = laplacian(Psi(:), Kx, Kz, Fr, alpha, beta, Lx, Lz);
% Xi2 = reshape(xi2, Nz, Nx);
combinedEigenfunctionPlot(X/Fr, Z, Xi, Psi, B)

%% 

% Xihat = fft2(Xi);
% Psihat = fft2(Psi);
% check_fftSymmetry(Xihat, Nx, Nz)
% check_fftSymmetry(Psihat, Nx, Nz)
% I = ones(Nz,Nx);
% kx = 2*pi/Lx;
% kz = 2*pi/Lz;
% L0 = -(Kx.^2 + Kz.^2);
% L  = -((Fr.^2).*((Kx + kx*alpha*I).^2) + (Kz + kz*beta*I).^2);
% Qhat        = zeros(size(Xihat));
% Qhat(L0~=0) = Xihat(L0~=0)./L(L0~=0);
% Qhat(L0==0) = 0;
% check_fftSymmetry(Qhat, Nx, Nz);
% Q           = real(ifft2(Qhat));
% combinedEigenfunctionPlot(X/Fr, Z, Xi, Q, B)