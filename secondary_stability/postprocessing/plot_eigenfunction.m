clear
close all
clc
 
%% ADD UTILITIES TO PATH

addpath('../../utility_belt');
addpath('..'); 

%% PLOT AESTHETICS

fs   = 16;
lw   = 3;
cmap = slanCM('magma');

%% PARAMETERS

Reb     = 100;
Fr      = 0.01;
alpha   = 0;
beta    = 0;
efn_idx = 6;

%% FILEPATH CONFIG - EIGENFUNCTION

ParentFolder = sprintf('../secondary_stability_solutions_maxit10000');
SolnFolder   = sprintf('Reb%.2f_alpha%.2f',Reb,alpha);
SolnFile     = sprintf('spectrum_Reb%.2f_alpha%.2f',Reb,alpha);
FilePath     = sprintf('%s/%s/%s',ParentFolder,SolnFolder,SolnFile);

%% LOAD REQUIRED EIGENSOLUTIONS DATA

load(FilePath);

%% GRAB EIGENVECTOR AND OTHER DATA

Nx = meta.resolution(1);
Nz = meta.resolution(2);
Lx = meta.domainsize(1);
Lz = meta.domainsize(2);
scale_factor = alpha;
if abs(scale_factor) == 0
    scale_factor = 1.0;
end
Lxp = Lx * scale_factor;
Lzp = Lz; 
kx  = (2*pi/Lxp) * [0:(Nx/2-1), (-Nx/2):-1];
kz  = (2*pi/Lzp) * [0:(Nz/2-1), (-Nz/2):-1];

[Kx, Kz] = meshgrid(kx, kz);

%% GENERATE 2D EIGENFUNCTION FIELDS (PLOT AND SAVE)

[Xip, Bp, Psip] = process_eigenfunction(eigvecs, efn_idx, Kx, Kz, Fr, alpha, beta, Lxp, Lzp);

%% PLOTTING BLOCK

xp   = (0:(Nx-1))*(Lxp/Nx);
zp   = (0:(Nz-1))*(Lzp/Nx);
efns = figure('WindowState', 'maximized', 'Color', 'w');
help_plot_3_fields(xp/Fr, zp, Xip, Bp, Psip, '$\xi^\prime$', '$b^\prime$', '$\psi^\prime$', cmap, fs, lw)
plotname = sprintf('Eigenfunctions_Reb%.2f_alpha%.2f_index%d.pdf',Reb,alpha,efn_idx);
exportgraphics(efns, plotname, 'ContentType', 'vector', 'Resolution', 500);
