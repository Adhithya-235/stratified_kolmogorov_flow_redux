function [eigvals, eigvecs, dom_mode, meta] = stability_eigensolve(Reb,...
    Pr, Fr, alpha, beta, numEig, residualtol, maxit, outputflag,...
    resolutionscale)

%% GET BASIC STATE

rescale = resolutionscale;
[Xihat, Bhat, Psihat, Lx, Lz, Nx, Nz] = nonlinear_bs(Reb, rescale);

%% START LOG

fprintf('Starting with Re_b = %f, Resolution = [%d, %d].\n', Reb, Nx, Nz)

%% FOURIER GRID

kxVals = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):-1];
kzVals = (2*pi/Lz)*[0:(Nz/2-1) (-Nz/2):-1];
[Kx, Kz] = meshgrid(kxVals, kzVals);

%% MATRIX OPERATOR AS A FUNCTION HANDLE

Afun = @(qhat) boussinesq(qhat, Xihat, Bhat, Psihat, Kx, Kz, alpha, beta,...
    Lx, Lz, Fr, Pr, Reb);

%% EIGENVALUE SOLVER OPTIONS

dim        = 2*Nx*Nz;
num        = numEig;    % COMPUTE 10 LARGEST EIGENVALUES (BY REAL PART)
type       = 'largestreal';
opts.tol   = residualtol;
opts.maxit = maxit;
opts.disp  = outputflag;
opts.p     = 51;
opts.v0    = construct_StartVector(Nx, Nz); 
opts.fail  = 'keep';

%% EIGENVALUE SOLVE

tIterStart = tic;
[V, D] = eigs(Afun, dim, num, type, opts);
tIter  = toc(tIterStart);

%% PROCESS DOMINANT EIGENFUNCTION

index = 1;
[Xi, B, Psi] = process_eigenfunction(V, index, Kx, Kz,  Fr, alpha, beta,...
    Lx, Lz);

%% PLOT DOMINANT EIGENFUNCTION

scale = alpha;
if scale == 0
    scale = 1;
end
Lxp      = (1/scale)*Lx;
Lzp      = Lz;
xp       = linspace(0, Lxp, Nx+1); xp = xp(1:end-1);
zp       = linspace(0, Lzp, Nz+1); zp = zp(1:end-1);

%% FUNCTION OUTPUTS

eigvals      = Fr*diag(D);
eigvecs      = V;
dom_mode.xp  = xp;
dom_mode.zp  = zp;
dom_mode.Xi  = Xi;
dom_mode.B   = B;
dom_mode.Psi = Psi;
meta.resolution = [Nx, Nz];
meta.domainsize = [Lxp, Lzp];
meta.itertime   = tIter; 

end