clear
close all
clc

%% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% FILE PARAMETERS

folder_name = '2024-09-29_02-47-09'; 
data_folder = 'results_branch4'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 7:9; 
wrap        = 0; 
unwrap      = 0; 

%% SIMULATION PARAMETERS

Fr = 0.02;
Rb = 50;
Pr = 1;
Lx = 1;
Lz = 4*pi/3;
Nx = 512;
Nz = 512;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DATA

[x, z, ~, ~]              = get_space_data(folder_name, data_folder, file_name, wrap);
[t, ~, ~, ~, vort, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%%  SELECT VORTICITY SLICE

tindex    = find(abs(t-8.49)<1e-3);
zindex    = find(abs(z-2.5)<dz);
vortslice = vort(zindex(1), :, tindex);

%% PLOT VORTICITY SLICE AS A FUNCTION OF x

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(x, vortslice, '-o', 'linewidth', 4, 'markersize', 3)
xlabel('$x$', 'interpreter', 'latex')
ylabel('$\xi$', 'interpreter', 'latex')
title(sprintf('Vorticity at t = %.2f and z = %.2f', t(tindex), z(zindex(1))))
set(gca, 'fontsize', 30)
axis tight
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
saveas(f1, 'vorticityslice.png') 

%% TAKE FFT

vorthat = fft(vortslice)/(Nx);
vortspc = vorthat.*conj(vorthat);

%% PLOT POWER SPECTRUM

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot((1:0.5*Nx)-1, vortspc(1:0.5*Nx), '-o', 'linewidth', 4, 'markersize', 5)
xlabel('Mode Number', 'interpreter', 'latex')
ylabel('$|\hat{\xi}|^2$', 'interpreter', 'latex')
title(sprintf('Enstrophy Spectrum at t = %.2f and z = %.2f', t(tindex), z(zindex(1))))
set(gca, 'fontsize', 30)
axis tight
xlim([0,49])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
saveas(f1, 'enstrophypspec.png') 
