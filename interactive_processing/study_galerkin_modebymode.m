clear
close all
clc

%% ADD PATH

addpath('../utility_belt'); 

%% FOLDER DETAILS

folder_name='2024-12-10_08-41-28'; 
data_folder='results_galerkin'; 
file_name='field_snapshots'; 

%% RUNTIME PARAMETERS

Nx=512; 
Nz=512; 
svec=1:50; 
modes=[0, 19, 23, 4]; 
wrap=0; 
unwrap=0; 
Fr=0.02; 

%% SIMULATION PARAMETERS

Rb = 50;
Pr = 1;
Lx = Fr*6*pi/0.34;
Lz = 4*pi/3;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DATA

[x, z, t, ~, ~, b, ~, nf] = get_galerkin_data(folder_name, data_folder, file_name, svec, Nx, Nz, modes, wrap);

%%  SELECT VORTICITY SLICE

zindex    = find(abs(z-2.5)<dz);
vortslice = squeeze(b(zindex(1), :, :));

%% TAKE FFT

vorthat = fftshift(fft(vortslice,[],1),1)/(Nx);
vortspc = vorthat.*conj(vorthat);

%% GET MODE NUMBERS

fmodes = fftshift([0:Nx/2, (-Nx/2+1):-1]);

%% PLOT LINE PLOT

f2 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, vortspc(257,:), '-o', 'linewidth', 3)
plot(t, vortspc(261,:), '-o', 'linewidth', 3)
plot(t, vortspc(276,:), '-o', 'linewidth', 3)
plot(t, vortspc(280,:), '-o', 'linewidth', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('PE Density', 'interpreter', 'latex')
legend('Mode 0','Mode 4','Mode 19', 'Mode 23', 'interpreter', 'latex', 'location', 'southeastoutside')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
ylim([1e-6, 1e-1])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'log')
drawnow
saveas(f2, 'lineplot-enstrophyspec-galerkin-2.png') 