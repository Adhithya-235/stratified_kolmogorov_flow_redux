clear
close all
clc

%% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% FILE PARAMETERS

folder_name = '2024-10-02_13-49-55'; 
data_folder = 'results_branch4'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 9; 
wrap        = 0; 
unwrap      = 0; 

%% SIMULATION PARAMETERS

Fr = 0.02;
Rb = 50;
Pr = 1;
Lx = Fr*6*pi/0.34;
Lz = 4*pi/3;
Nx = 1024;
Nz = 1024;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DATA

[x, z, X, Z]              = get_space_data(folder_name, data_folder, file_name, wrap);
[t, ~, ~, ~, vort, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%%  SELECT VORTICITY SLICE

tindex    = find(abs(t-8.45)<1e-3);
vortslice = vort(:, :, tindex);

%% PLOT VORTICITY AS A FUNCTION OF x AND z

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
subplot(211)
clim = [-3, 3];
help_plot_fields(X/Fr,Z,vortslice,clim,[],t(tindex))
title('Vorticity','interpreter','latex')

%% TAKE FFT, ZERO OUT MODE 19, RETURN

vorthat = fft(vortslice, [], 2)/(Nx);
vorthat(:,20) = zeros(Nz,1); vorthat(:,1006) = zeros(Nz,1);
vortmod = real(ifft(vorthat, [], 2));

%% PLOT MODIFIED VORTICITY SLICE

subplot(212)
clim = [min(min(vortmod)), max(max(vortmod))];
help_plot_fields(X/Fr,Z,vortmod,clim,[],t(tindex))
title('Vorticity Mode 19 Excised','interpreter','latex')
saveas(f1, 'vorticitysnapshot.png') 

