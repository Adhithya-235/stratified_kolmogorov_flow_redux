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
svec        = 21:50; 
wrap        = 0; 
unwrap      = 0; 

%% SIMULATION PARAMETERS

Fr = 0.02;
Rb = 50;
Pr = 1;
Lx = Fr*6*pi/0.34;
Lz = 4*pi/3;
Nx = 512;
Nz = 512;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DATA

[x, z, ~, ~]              = get_space_data(folder_name, data_folder, file_name, wrap);
[t, ~, ~, ~, vort, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% TAKE FFT

vorthat = fft(vort,[],2)/(Nx);
vortspc = vorthat.*conj(vorthat);

%% WRAP IN z AND TAKE AVERAGE

z         = [z; Lz];
vortspc   = [vortspc; vortspc(1,:,:)];
vortspct  = trapz(t, vortspc, 3)/(t(end)-t(1));
vortspczt = trapz(z, vortspct, 1)/Lz;

%% PLOT POWER SPECTRUM

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot((1:0.5*Nx)-1, vortspczt(1:0.5*Nx), '-o', 'linewidth', 4, 'markersize', 5)
xlabel('Mode Number', 'interpreter', 'latex')
ylabel('$|\hat{\xi}|^2$', 'interpreter', 'latex')
title(sprintf('Enstrophy Spectrum'))
set(gca, 'fontsize', 30)
axis tight
axis square
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'log', 'YScale', 'log')
saveas(f1, 'enstrophypspec.png') 

