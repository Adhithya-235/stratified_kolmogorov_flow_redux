%=======================================================================%
%   Plot hovmollers of the TKE budget. Best practice to run
%   from command line on a unix environment is to use the following
%   syntax:
%
%        matlab -batch "clear;close all;clc;addpath('../utility_belt'); ...
%               folder_name='$folder_name'; data_folder='$data_folder';...
%                    file_name='field_snapshots'; stride=20; svec=[1:3]; wrap=1; unwrap=0;...
%                          Ri = 0.1; Re = 10000; plot_tke_budget"
%
%=======================================================================%

%% READ DATA

[x, z, ~, ~]              = get_space_data(folder_name, data_folder, file_name, wrap);
[t, u, w, b, p, ~, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% GET DOMAIN LENGTH AND INTERVAL

Lx = x(end);
dx = x(2)-x(1);

%% CALCULATE MEAN AND PERTURBATION OF PRIMITIVE VARIABLES
disp('Calculating horizontal mean and anomaly of primitive variables.')
[um, up] = get_pert_fields(u,x,Lx,unwrap);
[~, wp]  = get_pert_fields(w,x,Lx,unwrap);
[~, bp]  = get_pert_fields(b,x,Lx,unwrap);
[~, pp]  = get_pert_fields(p,x,Lx,unwrap);

%% COMPUTE AUXILIARY FIELDS

disp('Calculating auxiliary fields.')
Ek          = 0.5*(up.^2 + wp.^2);
[dxu,dzu,~] = gradient(up,x,z,1);
[dxw,dzw,~] = gradient(wp,x,z,1);
dissfield   = (1/Re)*(dxu.^2 + dzu.^2 + dxw.^2 + dzw.^2);

%% COMPUTE AUXILIARY AVERAGES

disp('Calculating auxiliary averages.')
[tke,~]    = get_pert_fields(Ek,x,Lx,unwrap);
[wpbar,~]  = get_pert_fields(wp.*pp,x,Lx,unwrap);
[wEkbar,~] = get_pert_fields(wp.*Ek,x,Lx,unwrap);
[uwbar,~]  = get_pert_fields(wp.*up,x,Lx,unwrap);
[bwbar,~]  = get_pert_fields(wp.*bp,x,Lx,unwrap);

%% COMPUTE AUXILIARY GRADIENTS

disp('Calculating auxiliary gradients.')
[~,dzum,~] = gradient(um,1,z,1);
[~,dztk,~] = gradient(tke,1,z,1);

%% COMPUTE TKE BUDGET TERMS

disp('Computing TKE budget terms.')
[~,press_work,~] = gradient(-wpbar,1,z,1);  
[~,turb_trans,~] = gradient(-wEkbar,1,z,1);  
[~,shear_work,~] = gradient((1/Re)*dztk,1,z,1);   
shear_prod       = -uwbar.*dzum;
[tke_diss,~]     = get_pert_fields(-dissfield,x,Lx,unwrap);
buoy_prod        = Ri*bwbar;

%% PLOTTING PARAMETERS

yl = 4;

%% PLOT AND SAVE HOVMOLLERS -- TKE

disp('Plotting TKE hovmoller.')
f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
pcolor(t,z,squeeze(tke))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$z$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('Horz. Avg. TKE'))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
ylim([-yl,yl])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
drawnow
saveas(f1, sprintf('../%s/plots/tke_hovmollers/tke.fig', folder_name)) 
saveas(f1, sprintf('../%s/plots/tke_hovmollers/tke.png', folder_name)) 

%% PLOT AND SAVE HOVMOLLERS -- PRESSURE WORK

disp('Plotting pressure work hovmoller.')
f2 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
pcolor(t,z,squeeze(press_work))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$z$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('Horz. Avg. Pressure Work'))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
ylim([-yl,yl])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
drawnow
saveas(f2, sprintf('../%s/plots/tke_hovmollers/press_work.fig', folder_name)) 
saveas(f2, sprintf('../%s/plots/tke_hovmollers/press_work.png', folder_name)) 

%% PLOT AND SAVE HOVMOLLERS -- TURBULENT TRANSPORT

disp('Plotting turbulent transport hovmoller.')
f3 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
pcolor(t,z,squeeze(turb_trans))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$z$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('Horz. Avg. Turb. Transport'))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
ylim([-yl,yl])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
drawnow
saveas(f3, sprintf('../%s/plots/tke_hovmollers/turb_trans.fig', folder_name)) 
saveas(f3, sprintf('../%s/plots/tke_hovmollers/turb_trans.png', folder_name)) 

%% PLOT AND SAVE HOVMOLLERS -- SHEAR WORK

disp('Plotting shear work hovmoller.')
f4 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
pcolor(t,z,squeeze(shear_work))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$z$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('Horz. Avg. Shear Work'))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
ylim([-yl,yl])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
drawnow
saveas(f4, sprintf('../%s/plots/tke_hovmollers/shear_work.fig', folder_name)) 
saveas(f4, sprintf('../%s/plots/tke_hovmollers/shear_work.png', folder_name)) 

%% PLOT AND SAVE HOVMOLLERS -- SHEAR PRODUCTION

disp('Plotting shear production hovmoller.')
f5 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
pcolor(t,z,squeeze(shear_prod))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$z$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('Horz. Avg. Shear Production'))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
ylim([-yl,yl])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
drawnow
saveas(f5, sprintf('../%s/plots/tke_hovmollers/shear_prod.fig', folder_name)) 
saveas(f5, sprintf('../%s/plots/tke_hovmollers/shear_prod.png', folder_name)) 

%% PLOT AND SAVE HOVMOLLERS -- TKE DISSIPATION

disp('Plotting TKE dissipation hovmoller.')
f6 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
pcolor(t,z,squeeze(tke_diss))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$z$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('Horz. Avg. TKE Dissipation'))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
ylim([-yl,yl])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
drawnow
saveas(f6, sprintf('../%s/plots/tke_hovmollers/tke_diss.fig', folder_name)) 
saveas(f6, sprintf('../%s/plots/tke_hovmollers/tke_diss.png', folder_name)) 

%% PLOT AND SAVE HOVMOLLERS -- BUOYANT PRODUCTION

disp('Plotting buoyant production hovmoller.')
f7 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
pcolor(t,z,squeeze(buoy_prod))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$z$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('Horz. Avg. Buoyant Production'))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
ylim([-yl,yl])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)
drawnow
saveas(f7, sprintf('../%s/plots/tke_hovmollers/buoy_prod.fig', folder_name)) 
saveas(f7, sprintf('../%s/plots/tke_hovmollers/buoy_prod.png', folder_name)) 

%% PLOT AND SAVE -- RESOLUTION CHECK TIMESERIES

Lkmin    = ((Re^(-3))./(max(abs(squeeze(tke_diss))))).^(1/4);
rescheck = 2.5*Lkmin/dx;
disp('Plotting resolution check timeseries.')
f8 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
yline(1, 'g--', 'linewidth', 3)
plot(t, rescheck, '-o', 'linewidth', 3, 'markersize', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$2.5 L_K/\Delta x$', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
drawnow
saveas(f8, sprintf('../%s/plots/timeseries/resolution_check.fig', folder_name)) 
saveas(f8, sprintf('../%s/plots/timeseries/resolution_check.png', folder_name)) 

%% PLOT AND SAVE -- TKE RHS AND LHS

rhs       = press_work+turb_trans+shear_work+shear_prod+tke_diss+buoy_prod;
[~,~,lhs] = gradient(tke,1,1,t);

f9 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
pcolor(t,z,log10(abs(squeeze(lhs-rhs))))
colormap turbo
shading flat
c1 = colorbar;
ylabel('$z$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
title(sprintf('TKE Equation RHS'))
c1.FontSize = 30;
c1.Location = 'eastoutside';
axis tight
xlim([t(1), t(end)])
ylim([-yl,yl])
box on
set(gca, 'fontsize', 30, 'boxstyle', 'full', 'linewidth', 2)

saveas(f9, sprintf('../%s/plots/tke_hovmollers/rhs.fig', folder_name)) 
saveas(f9, sprintf('../%s/plots/tke_hovmollers/rhs.png', folder_name))