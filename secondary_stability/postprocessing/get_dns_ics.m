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
efn_idx = 1;

%% FILEPATH CONFIG - EIGENFUNCTION

ParentFolder = sprintf('../solutions_branch_1_Pr_1');
SolnFolder   = sprintf('Reb%.2f_alpha%.2f',Reb,alpha);
SolnFile     = sprintf('spectrum_Reb%.2f_alpha%.2f.mat',Reb,alpha);
FilePath     = sprintf('%s/%s/%s',ParentFolder,SolnFolder,SolnFile);

%% LOAD REQUIRED EIGENSOLUTIONS DATA

dataStruct = load(FilePath);

%% GRAB EIGENVECTOR AND OTHER DATA

V   = dataStruct.eigvecs;
Nx  = dataStruct.meta.resolution(1);
Nz  = dataStruct.meta.resolution(2);
Lxp = dataStruct.meta.domainsize(1); 
Lzp = dataStruct.meta.domainsize(2); 
kx  = (2*pi/Lxp) * [0:(Nx/2-1), (-Nx/2):-1];
kz  = (2*pi/Lzp) * [0:(Nz/2-1), (-Nz/2):-1];

[Kx, Kz] = meshgrid(kx, kz);

%% GENERATE 2D EIGENFUNCTION FIELDS (PLOT AND SAVE)

[Xip, Bp, Psip] = process_eigenfunction(V, efn_idx, Kx, Kz, Fr, alpha, beta, Lxp, Lzp);

%% INTERPOLATE FIELDS

Nx   = 256;
Nz   = 256;

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

x    = (0:(Nx-1))*(Lxp/Nx);
z    = (0:(Nz-1))*(Lzp/Nx);
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

H5FileName = sprintf('initialize_ecs_Reb%.2f_alpha%.2f_idx%d.h5',Reb,alpha,efn_idx);
datasets   = ["/zeta", "/psi", "/b"];
variables  = {Xi_in, Psi_in, B_in};

for idx = 1:length(datasets)
    h5create(H5FileName, datasets(idx), [Nx, Nz]);
    h5write(H5FileName, datasets(idx), variables{idx});
end

