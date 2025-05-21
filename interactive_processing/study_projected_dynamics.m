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
svec        = 1:15; 
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

%% GET FOURIER COEFFICIENTS OF vort, PROJECT ONTO 0, 4, 19, 23, IFFT

modes        = [0:Nx/2, (-Nx/2+1):-1];
vorth        = fft(vort,[],2);
projh        = zeros(size(vorth));

% MODE 0

projh(:,1,:) = vorth(:,1,:);

% MODE 4

projh(:,5,:)     = vorth(:,5,:);
projh(:,end-3,:) = vorth(:,end-3,:);

% MODE 19

projh(:,20,:)     = vorth(:,20,:);
projh(:,end-18,:) = vorth(:,end-18,:);

% MODE 23

projh(:,24,:)     = vorth(:,24,:);
projh(:,end-22,:) = vorth(:,end-22,:);

% IFFT

vortp = ifft(projh,[],2);

%% MAKE MOVIE AND COLLECT FRAMES

MOV = help_animate_fields(X/Fr,Z,vortp,'$\zeta_p(x,z,t)$',t,nf);

%% CREATE AVI FILE

fprintf('Writing movie to file...\n')
filename = sprintf('../%s/plots/movies/projected_vorticity_production.avi', folder_name);
v = VideoWriter(filename); 
v.FrameRate = 16;
open(v)
for i = 1:nf
   fprintf('Writing frame %d of %d.\n',i,nf)
   writeVideo(v, MOV(i)) 
end
close(v)

