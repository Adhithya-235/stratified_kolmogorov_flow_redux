%=======================================================================%
%   Animate snapshots for visualization purposes. Best practice to run
%   from command line on a unix environment is to use the following
%   syntax:
%
%        matlab -batch "clear;close all;clc; ...
%               folder_name='$folder_name'; svec=[1:3]; wrap=0; animate_fields"
%
%=======================================================================%

%% READ DATA

[x, z, t, ~, ~, ~, vort, nf] = get_galerkin_data(folder_name, data_folder, file_name, svec, nx, nz, modes, wrap);
[X, Z] = meshgrid(x', z);

%% MAKE MOVIE AND COLLECT FRAMES

MOV = help_animate_fields(X/Fr,Z,vort,'$\zeta(x,z,t)$',t,nf);

%% CREATE AVI FILE

fprintf('Writing movie to file...\n')
filename = sprintf('../%s/plots/movies/vorticity_production.avi', folder_name);
v = VideoWriter(filename); 
v.FrameRate = 16;
open(v)
for i = 1:nf
   fprintf('Writing frame %d of %d.\n',i,nf)
   writeVideo(v, MOV(i)) 
end
close(v)
