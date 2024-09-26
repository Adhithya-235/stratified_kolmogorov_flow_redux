%=======================================================================%
%   Plot volume averaged perturbation energy timeseries. Best practice to run
%   from command line on a unix environment is to use the following
%   syntax:
%
%        matlab -batch "clear;close all;clc;addpath('../utility_belt'); ...
%               folder_name='$folder_name'; data_folder='$data_folder';...
%                    file_name='field_snapshots'; stride=20; svec=[1:3]; wrap=1; unwrap=0;...
%                          Fr = 0.02; plot_pte"
%
%=======================================================================%

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

%tlin = 0:0.1:40;
%grcv = 0.00002*exp(2*0.07937*tlin);

%% PLOT TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, pke, '-o', 'linewidth', 3)
plot(t, ppe, '-o', 'linewidth', 3)
plot(t, pte, '-o', 'linewidth', 3)
%plot(tlin, grcv, 'ko', 'linewidth', 3)
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

saveas(f, sprintf('../%s/plots/timeseries/pertenergy_timeseries.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/pertenergy_timeseries.png', folder_name)) 
