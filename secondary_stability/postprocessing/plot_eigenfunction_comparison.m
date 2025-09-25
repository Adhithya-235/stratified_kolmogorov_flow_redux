clear
close all
clc
 
%% ADD UTILITIES TO PATH

addpath('../../utility_belt');
addpath('..'); 

%% PARAMETERS

Reb     = 100;
Fr      = 0.01;
alpha   = [0, 0.5, 0.33];
beta    = 0;
efn_idx = [6, 4, 2];
nalpha  = length(alpha); 

%% PLOT AESTHETICS

fs   = 16;
lw   = 3;
cmap = slanCM('magma');
figure('WindowState', 'maximized', 'Color', 'w')
efns = tiledlayout(nalpha,nalpha,'TileSpacing','compact','Padding','compact');
plotname = sprintf('eigenfunctions_Reb%.2f_oscillatory2_alter.png',Reb);
titlestr = ['Instability Modes (O2-ALT), ' '$Re_b = $' num2str(Reb)];
title(efns, titlestr, 'interpreter', 'latex');

%% PROCESSING AND PLOTTING LOOP

for idxalpha = 1:nalpha

    %% LOGGING

    fprintf('Working on alpha = %.2f.\n', alpha(idxalpha));

    %% FILEPATH CONFIG - EIGENFUNCTION

    ParentFolder = sprintf('../solutions_branch_1_Pr_1');
    SolnFolder   = sprintf('Reb%.2f_alpha%.2f',Reb,alpha(idxalpha));
    SolnFile     = sprintf('spectrum_Reb%.2f_alpha%.2f',Reb,alpha(idxalpha));
    FilePath     = sprintf('%s/%s/%s',ParentFolder,SolnFolder,SolnFile);

    %% LOAD REQUIRED EIGENSOLUTIONS DATA

    load(FilePath);

    %% GRAB EIGENVECTOR AND OTHER DATA

    Nx = meta.resolution(1);
    Nz = meta.resolution(2);
    Lx = meta.domainsize(1);
    Lz = meta.domainsize(2);
    kx = (2*pi/Lx) * [0:(Nx/2-1), (-Nx/2):-1];
    kz = (2*pi/Lz) * [0:(Nz/2-1), (-Nz/2):-1];

    [Kx, Kz] = meshgrid(kx, kz);

    %% GENERATE 2D EIGENFUNCTION FIELDS

    [Xip, ~, ~] = process_eigenfunction(eigvecs, efn_idx(idxalpha), Kx, Kz, Fr, alpha(idxalpha), beta, Lx, Lz);
    xp          = (0:(Nx-1))*(Lx/Nx);
    zp          = (0:(Nz-1))*(Lz/Nx);


    %% TILE VECTOR

    rownum  = 1 + (idxalpha - 1)*nalpha;
    tilevec = [1, idxalpha];

    %% PLOTTING 

    ax = nexttile(rownum, tilevec);
    pcolor(xp/Fr, zp, Xip); 
    colormap(cmap)
    shading interp;
    if idxalpha < nalpha
        ax.XTickLabel = [];
    else
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', fs-2);
    end
    ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs-2);
    ax.FontSize = fs - 4;
    ax.LineWidth = lw/2; 
    ax.TickLabelInterpreter = 'latex';
    daspect([1 1 1])

end

%% SAVE PLOT

exportgraphics(efns, plotname, 'ContentType', 'vector', 'Resolution', 500);