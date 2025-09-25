clear
close all
clc

% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% DNS FILE PARAMETERS

folder_name = '2025-09-04_08-02-34'; 
data_folder = 'results2_ecs100'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 1:15; 
wrap        = 0; 
unwrap      = 0; 

%% DNS SIMULATION PARAMETERS

Fr = 0.01;
Rb = 6;
Pr = 1;
Lx = 2*0.03;
Lz = 2*pi/3;
Nx = 128;
Nz = 128;
dx = Lx/Nx;
dz = Lz/Nz;
alpha = 0.5;

%% READ DNS DATA

[x, z, ~, ~]              = get_space_data(folder_name, data_folder, file_name, 1);
[t, ~, ~, b, vort, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% ECS FILE PARAMETERS

ecsFolder = sprintf('../secondary_stability/import');
ecsFile   = sprintf('ECS_real_field_for_Reb=%d.mat',Rb);
ecsPath   = sprintf('%s/%s',ecsFolder,ecsFile);

%% READ ECS DATA

ecs          = load(ecsPath);
ecs.omega    = interpft(interpft(repmat(ecs.omega, 1, 1/alpha), Nz, 1), Nx, 2);
ecs.Buoyancy = interpft(interpft(repmat(ecs.Buoyancy, 1, 1/alpha), Nz, 1), Nx, 2); 

%% PLOT ECS DATA

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
subplot(211)
clim = [-0.25, 0.25];
help_plot_fields(x(1:end-1)/Fr,z(1:end-1),ecs.Buoyancy,clim,[],t(1))
title('Buoyancy','interpreter','latex')
subplot(212)
clim = [-3, 3];
help_plot_fields(x(1:end-1)/Fr,z(1:end-1),ecs.omega,clim,[],t(1))
title('Vorticity','interpreter','latex')
saveas(f, 'initialguesssnapshot.png') 

%% CALCULATE PERTURBATIONS AND PERTURBATION ENERGIES

vortp = vort - ecs.omega;
bp    = b - ecs.Buoyancy;
enstr = vortp.^2;
poten = bp.^2;

%% PLOT PERTURBATION DATA

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
subplot(211)
clim = [-0.005, 0.005];
help_plot_fields(x(1:end-1)/Fr,z(1:end-1),bp(:,:,1),clim,[],t(1))
title('Buoyancy','interpreter','latex')
subplot(212)
clim = [-0.05, 0.05];
help_plot_fields(x(1:end-1)/Fr,z(1:end-1),vortp(:,:,1),clim,[],t(1))
title('Vorticity','interpreter','latex')
saveas(f, 'initialperturbnapshot.png') 

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

tlin1 = 0:0.001:1;
grcv1 = (gpoten(1))*exp(2*15*(tlin1-tlin1(1)));

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
% xlim([t(1), t(end)])
xlim([0, 1])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'log')
drawnow
    
%% SAVE PLOT

saveas(f, sprintf('../%s/plots/timeseries/pertenergy_timeseries_grver.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/pertenergy_timeseries_grver2.png', folder_name)) 
