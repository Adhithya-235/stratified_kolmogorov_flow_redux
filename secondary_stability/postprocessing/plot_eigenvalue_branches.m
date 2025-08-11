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

%% PREPROCESS TO FIND OK-REBS AND VALS

stat.Reb          = stat.Reb(~isnan(stat.tracked_eigs));
stat.tracked_eigs = stat.tracked_eigs(~isnan(stat.tracked_eigs));
osc1.Reb          = osc1.Reb(~isnan(osc1.tracked_eigs));
osc1.tracked_eigs = osc1.tracked_eigs(~isnan(osc1.tracked_eigs));
osc2.Reb          = osc2.Reb(~isnan(osc2.tracked_eigs));
osc2.tracked_eigs = osc2.tracked_eigs(~isnan(osc2.tracked_eigs));
osc3.Reb          = osc3.Reb(~isnan(osc3.tracked_eigs));
osc3.tracked_eigs = osc3.tracked_eigs(~isnan(osc3.tracked_eigs));

%% PLOT EIGENVALUES IN COMPLEX PLANE

trackedevals = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(real(stat.tracked_eigs), imag(stat.tracked_eigs), 'k-o', 'linewidth', lw)
plot(real(osc1.tracked_eigs), imag(osc1.tracked_eigs), 'b-o', 'linewidth', lw)
plot(real(osc1.tracked_eigs), -imag(osc1.tracked_eigs), 'b-o', 'linewidth', lw)
plot(real(osc2.tracked_eigs), imag(osc2.tracked_eigs), 'r-o', 'linewidth', lw)
plot(real(osc2.tracked_eigs), -imag(osc2.tracked_eigs), 'r-o', 'linewidth', lw)
plot(real(osc3.tracked_eigs), imag(osc3.tracked_eigs), 'm-o', 'linewidth', lw)
plot(real(osc3.tracked_eigs), -imag(osc3.tracked_eigs), 'm-o', 'linewidth', lw)
xlabel('$\mathrm{Fr}\,\mathcal{R}\left[\sigma\right]$', 'interpreter', 'latex')
ylabel('$\mathrm{Fr}\,\mathcal{I}\left[\sigma\right]$', 'interpreter', 'latex')
set(gca, 'fontsize', fs)
xlim([-0.1,0.15])
ylim([-1,1])
axis square
grid on
box on
set(gca, 'linewidth', lw, 'fontsize', fs, 'XScale', 'linear', 'YScale', 'linear', 'ticklabelinterpreter', 'latex', 'GridLineStyle', ':')
drawnow
plotname = sprintf('Alpha0.00_UnstableBranches.png');
exportgraphics(trackedevals, plotname, 'ContentType', 'vector', 'Resolution', 500);
