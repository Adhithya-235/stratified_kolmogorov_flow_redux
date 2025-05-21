clear
close all
clc

%% LOAD DATA

load('import/ECS_real_field_for_Reb=1.mat')

%% PLOT

[X, Z] = meshgrid(x, z);
fields = {Buoyancy, omega, psi};
titles = {'b', 'omega', 'psi'};
Fr     = 0.01;
figure;

for i = 1:3
    subplot(1,3,i);
    pcolor(X/Fr, Z, fields{i});
    shading interp;       % smooth shading
    axis equal tight;     % ensure tight fit and equal units
    daspect([1 1 1])
    colormap(parula);     % you can change this to your preference
    colorbar;
    title(titles{i});
end
