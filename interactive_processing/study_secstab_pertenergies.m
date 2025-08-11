clear
close all
clc

% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% DNS FILE PARAMETERS

folder_name = '2025-07-10_12-05-22'; 
data_folder = 'results2_ecs100'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 1:4; 
wrap        = 0; 
unwrap      = 0; 

%% DNS SIMULATION PARAMETERS

Fr = 0.01;
Rb = 100;
Pr = 1;
Lx = 0.03;
Lz = 2*pi/3;
Nx = 256;
Nz = 256;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DNS DATA

[x, z, ~, ~]              = get_space_data(folder_name, data_folder, file_name, 1);
[t, ~, ~, b, vort, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% ECS FILE PARAMETERS

ecsFolder = sprintf('../secondary_stability/import');
ecsFile   = sprintf('ECS_real_field_for_Reb=%d.mat',Rb);
ecsPath   = sprintf('%s/%s',ecsFolder,ecsFile);

%% READ ECS DATA

ecs          = load(ecsPath);
ecs.omega    = interpft(interpft(ecs.omega, Nz, 1), Nx, 2);
ecs.Buoyancy = interpft(interpft(ecs.Buoyancy, Nz, 1), Nx, 2); 

%% CALCULATE PERTURBATIONS AND PERTURBATION ENERGIES

vortp = vort - ecs.omega;
bp    = b - ecs.Buoyancy;
enstr = vortp.^2;
poten = bp.^2;

%% WRAP DATA BEFORE VOLUME AVERAGE

poten = cat(2, poten, poten(:, 1, :));
enstr = cat(2, enstr, enstr(:, 1, :));
poten = cat(1, poten, poten(1, :, :));
enstr = cat(1, enstr, enstr(1, :, :));
Lx    = x(end);
Lz    = z(end);

%% VOLUME AVERAGE PERTURBATION ENERGIES

disp('Starting volume average.')
genstr = calc_volm_avg(enstr,x,Lx,z,Lz);
disp('Done with enstrophy.')
gpoten = calc_volm_avg(poten,x,Lx,z,Lz);
disp('Ending volume average.')

%% CHECK AGAINST LINEAR GROWTH RATE

tlin1 = 0:0.005:2;
grcv1 = (5e-6)*exp(2*2.17*tlin1);

%% PLOT TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, genstr, '-o', 'linewidth', 3)
plot(t, gpoten, '-o', 'linewidth', 3)
plot(tlin1, grcv1, 'ko', 'linewidth', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('Pert Energy', 'interpreter', 'latex')
legend('Enstrophy','PE', 'SSGR', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'log')
drawnow
    
%% SAVE PLOT

saveas(f, sprintf('../%s/plots/timeseries/pertenergy_timeseries_grver.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/pertenergy_timeseries_grver.png', folder_name)) 
