clear
close all
clc

%% ADD UTILITY FUNCTIONS TO PATH

addpath('../utility_belt');

%% FILE PARAMETERS

folder_name = '2024-09-29_02-47-09'; 
data_folder = 'results_branch4'; 
file_name   = 'field_snapshots'; 
stride      = 1; 
svec        = 21:50; 
wrap        = 1; 
unwrap      = 0; 

%% SIMULATION PARAMETERS

Fr = 0.02;
Rb = 50;
Pr = 1;
Lx = Fr*6*pi/0.34;
Lz = 4*pi/3;
Nx = 512;
Nz = 512;
dx = Lx/Nx;
dz = Lz/Nz;

%% READ DATA

[x, z, ~, ~]           = get_space_data(folder_name, data_folder, file_name, wrap);
[t, u, ~, b, ~, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap);

%% GET FOURIER COEFFICIENTS OF vort, PROJECT ONTO 0, 4, 19, 23, IFFT

% modes        = [0:Nx/2, (-Nx/2+1):-1];
% uh           = fft(u,[],2);
% bh           = fft(b,[],2);
% projuh       = zeros(size(uh));
% projbh       = zeros(size(bh));

% % MODE 0

% projuh(:,1,:) = uh(:,1,:);
% projbh(:,1,:) = bh(:,1,:);

% % MODE 4

% projuh(:,5,:)     = uh(:,5,:);
% projuh(:,end-3,:) = uh(:,end-3,:);
% projbh(:,5,:)     = bh(:,5,:);
% projbh(:,end-3,:) = bh(:,end-3,:);

% % MODE 19

% projuh(:,20,:)     = uh(:,20,:);
% projuh(:,end-18,:) = uh(:,end-18,:);
% projbh(:,20,:)     = bh(:,20,:);
% projbh(:,end-18,:) = bh(:,end-18,:);

% % MODE 23

% projuh(:,24,:)     = uh(:,24,:);
% projuh(:,end-22,:) = uh(:,end-22,:);
% projbh(:,24,:)     = bh(:,24,:);
% projbh(:,end-22,:) = bh(:,end-22,:);

% % IFFT

% up = ifft(projuh,[],2);
% bp = ifft(projbh,[],2);

%% AVERAGE IN t AND x

ut  = trapz(t, up, 3)/(t(end)-t(1));
um  = trapz(x, ut, 2)/Lx;
bt  = trapz(t, bp, 3)/(t(end)-t(1));
bm  = trapz(x, bt, 2)/Lx;

%% PLOT DATA

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
subplot(131)
hold on
plot(um, z, '-o', 'linewidth', 4, 'markersize', 5)
xlabel('u mean', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
axis tight
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
subplot(132)
hold on
plot(z+bm, z, '-o', 'linewidth', 4, 'markersize', 5)
xlabel('b mean', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
axis tight
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'linear')
saveas(f1, 'means-fulldns.png') 

%% SAVE VARIABLES

save('means_projdns',"um","bm", "z");