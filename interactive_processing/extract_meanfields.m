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
wrap        = 1; 
unwrap      = 0; 

%% SIMULATION PARAMETERS

Fr = 0.02;
Rb = 50;
Pr = 1;
Lx = 1;
Lz = 4*pi/3;
Nx = 1024;
Nz = 1024;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DATA

[x, z, X, Z]           = get_space_data(folder_name, data_folder, file_name, wrap);
[t, u, ~, b, ~, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%%  SELECT VORTICITY SLICE

tindex = find(abs(t-8.14)<1e-3);
uslice = u(:, :, tindex);
bslice = b(:, :, tindex);

%% CALCULATE HORIZONTAL MEAN

[um, ~] = get_pert_fields(uslice,x,Lx,unwrap);
[bm, ~] = get_pert_fields(bslice,x,Lx,unwrap);

%% UNWRAP

um(end) = [];
bm(end) = [];
z(end)  = [];

%% INTERPOLATE 

Nz_target = 768;
dz_target = Lz/Nz_target;
z_target  = (0:Nz_target-1)*dz_target;
um_target = interpft(um,Nz_target);
bm_target = interpft(bm,Nz_target);

%% PLOT VARIABLES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
subplot(121)
hold on
plot(um_target, z_target, '-', 'linewidth', 4)
plot(sin(3*z_target), z_target, 'm--', 'linewidth', 4)
xlabel('$\overline{u}$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
axis tight
set(gca, 'fontsize', 30)
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
subplot(122)
plot(bm_target, z_target, '-', 'linewidth', 4)
xlabel('$\overline{b}$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
axis tight
set(gca, 'fontsize', 30)
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
saveas(f, 'linearstabsnapshot.png') 

%% SAVE DATA

save('linstab_basic',"um_target","bm_target","z_target","Rb","Fr","Pr","Lz","Nz_target","dz_target");
