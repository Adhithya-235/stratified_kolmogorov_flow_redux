clear
close all
clc

%% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% FILE PARAMETERS

folder_name = '2024-10-02_13-49-55'; 
data_folder = 'results_branch4'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 8:15; 
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

[x, z, X, Z]              = get_space_data(folder_name, data_folder, file_name, wrap);
[t, ~, ~, ~, vort, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%%  SELECT VORTICITY SLICE

zindex    = find(abs(z-2.5)<dz);
vortslice = squeeze(vort(zindex(1), :, :));

%% TAKE FFT

vorthat = fftshift(fft(vortslice,[],1),1)/(Nx);
vortspc = vorthat.*conj(vorthat);

%% GET MODE NUMBERS

modes = fftshift([0:Nx/2, (-Nx/2+1):-1]);

%% INITIALIZE FIGURE -- DENSITY PLOT

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% PLOT -- DENSITY PLOT

pcolor(t,modes,log10(vortspc))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$k$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('log10(Enstrophy), z = %.2f', z(zindex(1))))
caxis([-6 -0])
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
ylim([0,30])
xlim([t(1),t(end)])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
saveas(f, 'hovmoller-enstrophyspec.png') 

%% PLOT LINE PLOT

f2 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, vortspc(513,:), '-o', 'linewidth', 3)
plot(t, vortspc(517,:), '-o', 'linewidth', 3)
plot(t, vortspc(532,:), '-o', 'linewidth', 3)
plot(t, vortspc(536,:), '-o', 'linewidth', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('Enstrophy Density', 'interpreter', 'latex')
legend('Mode 0','Mode 4','Mode 19', 'Mode 23', 'interpreter', 'latex', 'location', 'southeastoutside')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
ylim([1e-6, 1e1])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'log')
drawnow
saveas(f2, 'lineplot-enstrophyspec.png') 