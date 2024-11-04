clear
close all
clc

%% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% FILE PARAMETERS

folder_name = '2024-10-18_15-01-21'; 
data_folder = 'results_branch1'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 1:15; 
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

[x, z, ~, ~]              = get_space_data(folder_name, data_folder, file_name, wrap);
[t, u, w, b, ~, nf]       = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% GET DOMAIN DIMENSIONS

Lx = x(end);
Lz = abs(z(end))+abs(z(1));

%% CALCULATE PERTURBATION FIELDS

[~, up] = get_pert_fields(u,x,Lx,unwrap);
[~, wp] = get_pert_fields(w,x,Lx,unwrap);
[~, bp] = get_pert_fields(b,x,Lx,unwrap);

%% FORM LOCAL KE AND TE

lke = up.^2 + (Fr*wp).^2;
lpe = (bp.^2);

%% CALCULATE VOLUME AVERAGE

disp('Starting volume average.')
pke = calc_volm_avg(lke,x,Lx,z,Lz);
disp('Done with kinetic energy.')
ppe = calc_volm_avg(lpe,x,Lx,z,Lz);
disp('Ending volume average.')
pte = pke + ppe;

%% INITIALIZE FIGURE

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% CHECK AGAINST LINEAR GROWTH RATE

tlin1 = 0:0.005:1;
grcv1 = 0.000005*exp(2*0.05*tlin1);
tlin2 = 1.75:0.005:2;
grcv2 = 0.000005*(exp(-2*20.5*1.75))*exp(2*20.5*tlin2);

%% PLOT TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, pke, '-o', 'linewidth', 3)
plot(t, ppe, '-o', 'linewidth', 3)
plot(t, pte, '-o', 'linewidth', 3)
plot(tlin1, grcv1, 'ko', 'linewidth', 3)
plot(tlin2, grcv2, 'mo', 'linewidth', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('Pert Energy', 'interpreter', 'latex')
legend('KE','PE','TE', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'log')
drawnow
    
%% SAVE PLOT

saveas(f, sprintf('../%s/plots/timeseries/pertenergy_timeseries_grver.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/pertenergy_timeseries_grver.png', folder_name)) 
