clear
close all
clc

%% SETUP TIMING START

tSetupStart = tic;

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

[Xihat, Bhat, Psihat] = laminar_bs(Z, m);

%% MATRIX OPERATOR AS A FUNCTION HANDLE

Afun = @(qhat) boussinesq(qhat, Xihat, Bhat, Psihat, Kx, Kz, alpha, beta,...
    Lx, Lz, Fr, Pr, Reb);
% targeteig = (0.2561 + 1i*0)/Fr;
% Afun = @(qhat) linop(qhat)-targeteig*qhat;

%% EIGENVALUE SOLVER OPTIONS

dim        = 2*Nx*Nz;
num        = 1;    % COMPUTE 10 LARGEST EIGENVALUES (BY REAL PART)
type       = 'largestreal';
opts.tol   = 1e-12;
opts.maxit = 2000;
opts.disp  = 1;
opts.p     = 20;
opts.v0    = construct_StartVector(Nx, Nz); 
opts.fail  = 'keep';

%% SETUP TIMER STOP

tSetup = toc(tSetupStart);

%% EIGENVALUE SOLVE

tIterStart = tic;
[V, D] = eigs(Afun, dim, num, type, opts);
tIter  = toc(tIterStart);
disp('Computed eigenvalues:');
disp(Fr*diag(D))

%% PROCESS AND PLOT DOMINANT EIGENFUNCTION

index = 1;
[Xi, B, Psi] = process_eigenfunction(V, index, Kx, Kz,  Fr, alpha, beta,...
    Lx, Lz);
combinedEigenfunctionPlot(X/Fr, Z, Xi, Psi, B)

%% PRINT TIMES

fprintf('Setup time: %.4f seconds\n', tSetup);
fprintf('Iteration time: %.4f seconds\n', tIter);
