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
svec=21:50; 
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

[x, z, t, u, ~, b, ~, nf] = get_galerkin_data(folder_name, data_folder, file_name, svec, Nx, Nz, modes, wrap);

%% AVERAGE IN t AND x

ut  = trapz(t, u, 3)/(t(end)-t(1));
um  = trapz(x, ut, 2)/Lx;
bt  = trapz(t, b, 3)/(t(end)-t(1));
bm  = trapz(x, bt, 2)/Lx;

%% PLOT DATA

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
subplot(131)
hold on
plot(um, z, '-o', 'linewidth', 4, 'markersize', 5)
xlabel('u mean', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
axis tight
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
subplot(132)
hold on
plot(z+bm, z, '-o', 'linewidth', 4, 'markersize', 5)
xlabel('b mean', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
axis tight
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
saveas(f1, 'means-fulldns.png') 

%% SAVE VARIABLES

save('means_galerk',"um","bm", "z");