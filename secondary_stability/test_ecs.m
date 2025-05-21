clear
close all
clc

%% PHYSICAL PARAMETERS

Reb   = 1;                  % Buoyancy Reynolds
Pr    = 1;                  % Prandtl
Fr    = 0.01;               % Froude
alpha = 0;                  % x Floquet Modifier
beta  = 0;                  % z Floquet Modifier

%% SELECT OPTIMAL SCALE FOR GRID CONVERGENCE

residualtol = 1e-6;
maxit       = 2000;
optscale = select_optimal_grid(Reb, Pr, Fr, alpha, beta, residualtol,...
    maxit);

%% SOLVE FULL PROBLEM

numEig      = 1;
residualtol = 1e-6;
maxit       = 2000;
outputflag  = 1;
[eigvals, eigvecs, dom_mode, meta] = stability_eigensolve(Reb, Pr, Fr,...
    alpha, beta, numEig, residualtol, maxit, outputflag, optscale);

%% LIST EIGENVALUES

disp('Computed eigenvalues:');
disp(eigvals);

%% PLOT DOMINANT EIGENFUNCTION

[Xp, Zp] = meshgrid(dom_mode.xp, dom_mode.zp);

combinedEigenfunctionPlot(Xp/Fr, Zp, dom_mode.Xi, dom_mode.Psi, dom_mode.B)

%% PRINT TIMES

fprintf('Iteration time: %.4f seconds\n', meta.itertime);
