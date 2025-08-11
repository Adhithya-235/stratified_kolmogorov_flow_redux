clear
close all
clc

%% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% FILE PARAMETERS

folder_name = '2025-07-07_06-28-48'; 
data_folder = 'results_ecs100'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 11:30; 
wrap        = 1; 
unwrap      = 0; 

%% SIMULATION PARAMETERS

Fr = 0.01;
Rb = 100;
Pr = 1;
Lx = 0.03;
Lz = 2*pi/3;
Nx = 256;
Nz = 256;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DATA

[x, z, ~, ~]           = get_space_data(folder_name, data_folder, file_name, wrap);
[t, ~, ~, ~, vort, nf] = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% FORM LOCAL ENSTROPHY

lens = (vort.^2)/2;

%% CALCULATE VOLUME AVERAGE

disp('Starting DNS volume averages.')
pens = calc_volm_avg(lens,x,Lx,z,Lz);
disp('Done with enstrophy. Ending DNS volume averages.')

%% ECS FILE PARAMETERS

ecsFolder = sprintf('../secondary_stability/import');
ecsFile   = sprintf('ECS_real_field_for_Reb=%d.mat',Rb);
ecsPath   = sprintf('%s/%s',ecsFolder,ecsFile);

%% READ ECS DATA

ecs    = load(ecsPath);
ecsV   = ecs.omega;
ecsV   = interpft(interpft(ecsV, Nz, 1), Nx, 2);

%% WRAP ECS AND COMPUTE POTENTIAL ENERGY

ecsV   = cat(2, ecsV, ecsV(:, 1, :));
ecsV   = cat(1, ecsV, ecsV(1, :, :));  
ecsens = (ecsV.^2)/2;
disp('Computing ECS Enstrophy.')
ecsEns = calc_volm_avg(ecsens,x,Lx,z,Lz);
disp('Finished')

%% PLOT ENSTROPHY COMPARISON TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, pens, '-o', 'linewidth', 3)
yline(ecsEns, 'k--', 'linewidth', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('Enstrophy', 'interpreter', 'latex')
legend('DNS','ECS', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
drawnow
saveas(f, sprintf('../%s/plots/timeseries/enstrophy_timeseries.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/enstrophy_timeseries.png', folder_name)) 

