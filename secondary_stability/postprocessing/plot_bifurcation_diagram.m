clear
close all
clc

%% PLOT AESTHETICS

fs   = 24;
lw   = 3;

%% LOAD FILES

stat = load('stationary_mode_data.mat');
osc1 = load('oscillatory_mode_1_data.mat');
osc2 = load('oscillatory_mode_2_data.mat');
osc3 = load('oscillatory_mode_3_data.mat');
Reb  = stat.Reb;

%% PREPROCESS TO FIND OK-REBS AND VALS

stat.tracked_eigs(real(stat.tracked_eigs)<0) = NaN;
osc1.tracked_eigs(real(osc1.tracked_eigs)<0) = NaN;
osc2.tracked_eigs(real(osc2.tracked_eigs)<0) = NaN;
osc3.tracked_eigs(real(osc3.tracked_eigs)<0) = NaN;

%% PLOT EIGENVALUES IN COMPLEX PLANE

trackedevals = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(stat.Reb, real(stat.tracked_eigs), '-o', 'linewidth', lw)
plot(osc1.Reb, real(osc1.tracked_eigs), '-s', 'linewidth', lw)
plot(osc2.Reb, real(osc2.tracked_eigs), '-s', 'linewidth', lw)
plot(osc3.Reb, real(osc3.tracked_eigs), '-s', 'linewidth', lw)
xlabel('$\mathrm{Re}_b$', 'interpreter', 'latex')
ylabel('$\mathrm{Fr}\,\mathcal{R}\left[\sigma\right]$', 'interpreter', 'latex')
legend('stat', 'osc1', 'osc2', 'osc3')
set(gca, 'fontsize', fs)
xlim([Reb(1), Reb(end)])
axis square
grid on
box on
set(gca, 'linewidth', lw, 'fontsize', fs, 'XScale', 'log', 'YScale', 'linear', 'ticklabelinterpreter', 'latex', 'GridLineStyle', ':')
drawnow
plotname = sprintf('Alpha0.50_Instabilities.png');
exportgraphics(trackedevals, plotname, 'ContentType', 'vector', 'Resolution', 500);
