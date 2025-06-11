clear
close all
clc

%% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% FILE PARAMETERS

folder_name = '2025-06-10_15-19-27'; 
data_folder = 'results_ecs'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 20; 
wrap        = 0; 
unwrap      = 0; 

%% SIMULATION PARAMETERS

Fr = 0.01;
Rb = 10;
Pr = 1;
Lx = 0.03;
Lz = 2*pi/3;
Nx = 128;
Nz = 128;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DATA

[x, z, X, Z]              = get_space_data(folder_name, data_folder, file_name, wrap);
[t, ~, ~, b, vort, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%%  SELECT VORTICITY SLICE

tindex    = find(abs(t-19.99)<1e-3);
vortslice = vort(:, :, tindex);
bslice    = b(:, :, tindex);

%% PLOT VARIABLES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
subplot(211)
clim = [0, 2*pi/3];
help_plot_fields(X/Fr,Z,z+bslice,clim,[],t(tindex))
title('Buoyancy','interpreter','latex')
subplot(212)
clim = [-3, 3];
help_plot_fields(X/Fr,Z,vortslice,clim,[],t(tindex))
title('Vorticity','interpreter','latex')
saveas(f, 'initialguesssnapshot.png') 

%% SAVE DATA

% save('initial_guess_adhi',"vortslice","bslice", "x","z","Rb","Fr","Pr","Lx","Lz","Nx","Nz","dx","dz");
