%=======================================================================%
%   Plot field snapshots for visualization purposes. Best practice to run
%   from command line on a unix environment is to use the following
%   syntax:
%
%        matlab -batch "clear;close all;clc; ...
%               folder_name='$folder_name'; svec=[1:3]; wrap=0; plot_fields"
%
%=======================================================================%

%% READ DATA

[x, z, X, Z]           = get_space_data(folder_name, data_folder, file_name, wrap);
[t, u, ~, b, vort, nf] = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% INITIALIZE FIGURE

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% MAKE AND SAVE PLOT

for i = 1:nf

    %% PLOT b

    subplot(311)
    clim = [0, 4*pi/3];
    help_plot_fields(X/Fr,Z,z+b(:,:,i),clim,[],t(i))
    title('Buoyancy','interpreter','latex')

    %% PLOT u

    subplot(312)
    clim = [-1, 1];
    help_plot_fields(X/Fr,Z,u(:,:,i),clim,[],t(i))
    title('Horz. Velocity','interpreter','latex')

    %% PLOT vort

    subplot(313)
    %clim = [min(min(min(vort))), max(max(max(vort)))];
    clim = [-3, 3];
    help_plot_fields(X/Fr,Z,vort(:,:,i),clim,[],t(i))
    title('Vorticity','interpreter','latex')

    %% DRAW FRAME

    drawnow
    
    %% SAVE PLOT
     
    saveas(f, sprintf('../%s/plots/frames/frame_time_%3.3f.png', folder_name, t(i))) 

end