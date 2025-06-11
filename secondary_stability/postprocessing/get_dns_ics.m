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

Reb     = 13;
Fr      = 0.01;
alpha   = 0;
beta    = 0;
efn_idx = 1;

%% FILEPATH CONFIG - EIGENFUNCTION

ParentFolder = sprintf('../secondary_stability_solutions_maxit10000');
SolnFolder   = sprintf('Reb%.2f_alpha%.2f',Reb,alpha);
SolnFile     = sprintf('converged_Reb%dp00_alpha%dp00.mat',Reb,alpha);
VarName      = sprintf('converged_Reb%dp00_alpha%dp00',Reb,alpha);
FilePath     = sprintf('%s/%s/%s',ParentFolder,SolnFolder,SolnFile);

%% LOAD REQUIRED EIGENSOLUTIONS DATA

eigenSoln = load(FilePath);
dataStruct = eigenSoln.(VarName);

%% GRAB EIGENVECTOR AND OTHER DATA

V  = dataStruct.converged_eigenvectors;
Nx = dataStruct.original_resolution(1);
Nz = dataStruct.original_resolution(2);
Lx = dataStruct.original_domainsize(1); 
Lz = dataStruct.original_domainsize(2); 
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

[Xip, Bp, Psip] = process_eigenfunction(V, efn_idx, Kx, Kz, Fr, alpha, beta, Lxp, Lzp);

%% INTERPOLATE FIELDS

Nx   = 128;
Nz   = 128;

Xip  = interpft(interpft(Xip, Nz, 1), Nx, 2);
Psip = interpft(interpft(Psip, Nz, 1), Nx, 2);
Bp   = interpft(interpft(Bp, Nz, 1), Nx, 2); 

%% NORMALIZE 

normfac = 1e3;
Xip     = normfac*Xip;
Bp      = normfac*Bp;
Psip    = normfac*Psip;

%% LOAD ECS

[Xi, B, Psi] = get_ecs_fields(Reb, Nx/128);

%% CONSTRUCT INITIAL CONDITION

Xi_in  = Xi + Xip;
B_in   = B + Bp;
Psi_in = Psi + Psip;

%% PLOTTING BLOCK -- OPTIONAL, ONLY FOR DEBUGGINH

x    = (0:(Nx-1))*(Lx/Nx);
z    = (0:(Nz-1))*(Lz/Nx);
xp   = (0:(Nx-1))*(Lxp/Nx);
zp   = (0:(Nz-1))*(Lzp/Nx);

efns = figure('WindowState', 'maximized', 'Color', 'w');
help_plot_3_fields(xp/Fr, zp, Xip, Bp, Psip, '$\xi^\prime$', '$b^\prime$', '$\psi^\prime$', cmap, fs, lw)
exportgraphics(efns, 'efns.pdf', 'ContentType', 'vector', 'Resolution', 500);

ecs = figure('WindowState', 'maximized', 'Color', 'w');
help_plot_3_fields(x/Fr, z, Xi, B, Psi, '$\xi_s$', '$b_s$', '$\psi_s$', cmap, fs, lw)
exportgraphics(ecs, 'ecs.pdf', 'ContentType', 'vector', 'Resolution', 500);

init = figure('WindowState', 'maximized', 'Color', 'w');
help_plot_3_fields(xp/Fr, zp, Xi_in, B_in, Psi_in, '$\xi_{in}$', '$b_{in}$', '$\psi_{in}$', cmap, fs, lw)
exportgraphics(init, 'init.pdf', 'ContentType', 'vector', 'Resolution', 500);

%% SAVE AS HDF5

H5FileName = sprintf('initialize_ecs_Reb%.2f_alpha%.2f.h5',Reb,alpha);
datasets   = ["/zeta", "/psi", "/b"];
variables  = {Xi_in, Psi_in, B_in};

for idx = 1:length(datasets)
    h5create(H5FileName, datasets(idx), [Nx, Nz]);
    h5write(H5FileName, datasets(idx), variables{idx});
end

