clear
close all
clc

%% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% FILE PARAMETERS

folder_name = '2025-08-19_17-48-03'; 
data_folder = 'results2_ecs100'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 1:25; 
wrap        = 1; 
unwrap      = 0; 

%% SIMULATION PARAMETERS

Fr = 0.01;
Rb = 6;
Pr = 1;
Lx = 2*0.03;
Lz = 2*pi/3;
Nx = 128;
Nz = 128;
dx = Lx/Nx;
dz = Lz/Nz;
alpha = 0.05;

%% READ DATA


[x, z, ~, ~]           = get_space_data(folder_name, data_folder, file_name, wrap);
[t, u, w, b, vort, nf] = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% FORM LOCAL KE AND TE

lke  = (u.^2 + (Fr*w).^2)/2;
lpe  = (b.^2)/2;
lens = (vort.^2)/2;

%% CALCULATE VOLUME AVERAGE

disp('Starting DNS volume averages.')
pke  = calc_volm_avg(lke,x,Lx,z,Lz);
disp('Done with kinetic energy.')
ppe  = calc_volm_avg(lpe,x,Lx,z,Lz);
disp('Done with potential energy.')
pens = calc_volm_avg(lens,x,Lx,z,Lz);
disp('Done with enstrophy. Ending DNS volume averages.')
pte  = pke + ppe;

%% INITIALIZE FIGURE

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% PLOT TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, pke, '-o', 'linewidth', 3)
plot(t, ppe, '-o', 'linewidth', 3)
plot(t, pte, '-o', 'linewidth', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('Energy', 'interpreter', 'latex')
legend('KE','PE','TE', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
drawnow
    
%% SAVE TIMESERIES PLOT

saveas(f, sprintf('../%s/plots/timeseries/energy_timeseries.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/energy_timeseries.png', folder_name)) 

%% ECS FILE PARAMETERS

ecsFolder = sprintf('../secondary_stability/import');
ecsFile   = sprintf('ECS_real_field_for_Reb=%d.mat',Rb);
ecsPath   = sprintf('%s/%s',ecsFolder,ecsFile);

%% READ ECS DATA

ecs    = load(ecsPath);
ecsB   = interpft(interpft(repmat(ecs.Buoyancy, 1, 1/alpha), Nz, 1), Nx, 2); 
ecsV   = interpft(interpft(repmat(ecs.omega, 1, 1/alpha), Nz, 1), Nx, 2);
ecsPsi = interpft(interpft(repmat(ecs.psi, 1, 1/alpha), Nz, 1), Nx, 2); 

%% COMPUTE u AND w FROM Psi

kx       = (2*pi/Lx)*[0:Nx/2, (-Nx/2+1):-1];
kz       = (2*pi/Lz)*[0:Nz/2, (-Nz/2+1):-1];
[Kx, Kz] = meshgrid(kx,kz);
ecsU     = real(ifft2(1i*Kz.*fft2(ecsPsi)));
ecsW     = -real(ifft2(1i*Kx.*fft2(ecsPsi)));

%% WRAP ECS AND COMPUTE POTENTIAL ENERGY

ecsV   = cat(2, ecsV, ecsV(:, 1, :));
ecsV   = cat(1, ecsV, ecsV(1, :, :));  
ecsU   = cat(2, ecsU, ecsU(:, 1, :));
ecsU   = cat(1, ecsU, ecsU(1, :, :));
ecsW   = cat(2, ecsW, ecsW(:, 1, :));
ecsW   = cat(1, ecsW, ecsW(1, :, :));
ecsB   = cat(2, ecsB, ecsB(:, 1, :));
ecsB   = cat(1, ecsB, ecsB(1, :, :));
ecste  = (ecsU.^2 + (Fr*ecsW).^2 + ecsB.^2)/2;
ecsens = (ecsV.^2)/2;
disp('Computing ECS TE.')
ecsTE = calc_volm_avg(ecste,x,Lx,z,Lz);
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

%% PLOT TE COMPARISON TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, pte, '-o', 'linewidth', 3)
yline(ecsTE, 'k--', 'linewidth', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('Total Energy', 'interpreter', 'latex')
legend('DNS','ECS', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
drawnow
saveas(f, sprintf('../%s/plots/timeseries/total_energy_timeseries.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/total_energy_timeseries.png', folder_name)) 
