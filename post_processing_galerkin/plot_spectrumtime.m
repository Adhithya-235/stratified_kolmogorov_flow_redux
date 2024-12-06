%=======================================================================%
%   Plot interfacial buoyancy spectra as a function of time. 
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
Nx = length(x);

%% CALCULATE PERTURBATION BUOYANCY

[~, bp] = get_pert_fields(b,x,Lx,unwrap);

%% GET INTERFACE INDICES

top_index = find(abs(z-2)<=1e-1);
top_index = top_index(end);

%% SLICE PERTURBATION BUOYANCY

bp_top = squeeze(bp(top_index,:,:));

%% CALCULATE POWER SPECTRUM

bptop_coeffs = fftshift((1/Nx)*fft(bp_top,[],1),1);

%% CALCULATE POWER

bptop_pow = bptop_coeffs.*conj(bptop_coeffs);

%% GET MODE NUMBERS

modes = fftshift([0:Nx/2, (-Nx/2+1):-1]);

%% INITIALIZE FIGURE -- DENSITY PLOT

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% PLOT -- DENSITY PLOT

pcolor(t,modes,log10(bptop_pow))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$\alpha$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('log10(PE), z = %.2f', z(top_index)))
caxis([-3 0])
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
ylim([0,10])
xlim([t(1),t(end)])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
daspect([6,1,1])
drawnow

%% SAVE PLOT -- DENSITY PLOT

saveas(f, sprintf('../%s/plots/timeseries/interfacial_pespectrum_timeseries.fig', folder_name)) 
saveas(f, sprintf('../%s/plots/timeseries/interfacial_pespectrum_timeseries.png', folder_name)) 

%% INITIALIZE FIGURE -- LINE PLOT

% f2 = figure;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
% hold on

%% PLOT -- LINE PLOT

% plot(modes,bptop_pow(:,1:10:100),'b-o','linewidth',4)%
% xlabel('$\alpha$', 'interpreter', 'latex')
% ylabel('PE Spectrum', 'interpreter', 'latex')
% title(sprintf('PE, z = %.2f', z(top_index)))
% axis tight
% xlim([0 modes(end)])
% ylim([1e-6 1e2])
% grid on
% box on
% set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2,'xscale','log','yscale','log')
% drawnow

%% SAVE PLOT -- LINE PLOT

% saveas(f2, sprintf('../%s/plots/timeseries/interfacial_pespectrum_lineplot.png', folder_name)) 