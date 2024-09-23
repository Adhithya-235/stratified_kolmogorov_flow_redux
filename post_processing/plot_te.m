%=======================================================================%
%   Plot a timeseries of cross-wind kinetic energy. Best practice to run
%   from command line on a unix environment is to use the following
%   syntax:
%
%        matlab -batch "clear;close all;clc; ...
%               folder_name='$folder_name'; maxs=1; plot_cwke"
%
%=======================================================================%


%% READ DATA

[t, ~, ke, pe, te] = get_timeseries_data(data_folder, folder_name, maxs);

%% PLOT TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, ke, '-', 'linewidth', 4)
plot(t, pe, '-', 'linewidth', 4)
plot(t, te, '-', 'linewidth', 4)
xlabel('$t$', 'interpreter', 'latex')
ylabel('Energy', 'interpreter', 'latex')
legend('KE','PE','TE', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')

%% SAVE VARIABLES AND PLOT

saveas(f, sprintf('../%s/plots/timeseries/te_timeseries.png', folder_name)) 

