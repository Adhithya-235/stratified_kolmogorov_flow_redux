%=======================================================================%
%   Plot interfacial buoyancy anomaly as a function of time. 
%   Best practice to run from command line on a unix environment 
%   is to use the following syntax:
%
%        matlab -batch "clear;close all;clc;addpath('../utility_belt'); ...
%               folder_name='$folder_name'; data_folder='$data_folder';...
%                    file_name='field_snapshots'; stride=20; svec=[1:3]; wrap=1;...
%                    unwrap=1; Ri = 0.1; plot_spectrumtime"
%
%=======================================================================%

%% READ DATA

[x, z, ~, ~]              = get_space_data(folder_name, data_folder, file_name, wrap);
[t, ~, ~, b, ~, ~, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% GET DOMAIN DIMENSIONS

Lx = x(end);

%% CALCULATE PERTURBATION BUOYANCY

[~, bp] = get_pert_fields(b,x,Lx,unwrap);

%% GET INTERFACE INDICES

top_index = find(abs(z-1.5)<=1e-1);
top_index = top_index(end);

%% SLICE PERTURBATION BUOYANCY

bp_top = squeeze(bp(top_index,:,:));

%% INITIALIZE FIGURE -- DENSITY PLOT

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% PLOT -- DENSITY PLOT

pcolor(t,x,bp_top)
colormap turbo
shading flat
c1 = colorbar;
ylabel('$x$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('Pert. b, z = %.2f', z(top_index)))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
drawnow

%% SAVE PLOT -- DENSITY PLOT

saveas(f, sprintf('../%s/plots/timeseries/interfacial_buoyancy_hovmoller.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/interfacial_buoyancy_hovmoller.png', folder_name)) 

