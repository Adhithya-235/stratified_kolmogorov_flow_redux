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
alpha   = 0.5;
beta    = 0;
efn_idx = 2;

%% FILEPATH CONFIG - EIGENFUNCTION

ParentFolder = sprintf('../solutions_branch_1_Pr_1');
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
Lxs = Lx * alpha;
if alpha == 0
    Lxs = Lx;
end
Lzs = Lz; 
kx  = (2*pi/Lxs) * [0:(Nx/2-1), (-Nx/2):-1];
kz  = (2*pi/Lzs) * [0:(Nz/2-1), (-Nz/2):-1];

[Kx, Kz] = meshgrid(kx, kz);

%% GENERATE 2D EIGENFUNCTION FIELDS (PLOT AND SAVE)

[Xip, Bp, Psip] = process_eigenfunction_2(eigvecs, efn_idx, Kx, Kz, Fr, alpha, beta, Lxs, Lzs);

%% PLOTTING BLOCK

xp   = (0:(Nx-1))*(Lx/Nx);
zp   = (0:(Nz-1))*(Lz/Nx);
efns = figure('WindowState', 'maximized', 'Color', 'w');
help_plot_3_fields(xp/Fr, zp, Xip, Bp, Psip, '$\xi^\prime$', '$b^\prime$', '$\psi^\prime$', cmap, fs, lw)
plotname = sprintf('Eigenfunctions_Reb%.2f_alpha%.2f_index%d.png',Reb,alpha,efn_idx);
exportgraphics(efns, plotname, 'ContentType', 'vector', 'Resolution', 500);
