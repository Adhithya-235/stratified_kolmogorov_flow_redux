clear 
close all
clc

%% REB VECTOR

Reb = [1,2,10,11,12,13,15,20,30,50,55,60,85,90,100,130,150,170,185,200];

%% PREALLOCATE FRAMES

MOV(length(Reb)) = struct('cdata',[],'colormap',[]);

%% PLOT DETAILS

fs = 30;
lw = 3;
f = figure('WindowState','maximized');

%% MAIN LOOP

for i = 1:length(Reb)

%% LOAD DATA

Reb(i)
load(sprintf('ECS_real_field_for_Reb=%d.mat', Reb(i)))

%% GET Nx, Nz, Lx, Lz

nx = length(x); nz = length(z);
Lx = 2*x(end) - x(end-1);
Lz = 2*z(end) - z(end-1);

%% PLOT OMEGA IN REAL SPACE 

subplot(121)
pcolor(x/0.01,z,omega)
shading interp;
axis tight;
daspect([1 1 1])
colormap turbo;
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.FontSize = fs;
c.Label.String = '$\omega$';
xlabel('$x/Fr$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs);
title(sprintf('ECS Vorticity, Reb =  %d', Reb(i)), 'FontSize', fs);
box on;
set(gca, 'LineWidth', lw, 'FontSize', fs, 'BoxStyle', 'full');

%% 2D FFT

omegahat  = (1/(nx*nz))*fftshift(fftshift(fft2(omega),1),2);
enstrophy = omegahat.*conj(omegahat)/2;
enstrophy = enstrophy./(max(max(enstrophy)));

%% WAVENUMBER GRID

kx = [0:nx/2, (-nx/2+1):-1];
kz = [0:nz/2, (-nz/2+1):-1];
kx = fftshift(kx);
kz = fftshift(kz);
[KX, KZ] = meshgrid(kx, kz);

%% PLOT ENSTROPHY

subplot(122)
pcolor(KX,KZ,log10(enstrophy))
shading faceted;
axis tight square;
colormap turbo;
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.FontSize = fs;
c.Label.String = 'log10($\hat{\omega}^2/2$)';
xlim([-40, 40])
ylim([-40, 40])
clim([-6, 1])
xlabel('$n_x$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$n_z$', 'Interpreter', 'latex', 'FontSize', fs);
title('ECS Enstrophy', 'Interpreter', 'latex', 'FontSize', fs);
box on;
set(gca, 'LineWidth', lw, 'FontSize', fs, 'BoxStyle', 'full');
drawnow
MOV(i) = getframe(gcf);

end

%% SAVE MOVIE

filename = 'enstrophy.avi';
v = VideoWriter(filename);
v.FrameRate = 12;
open(v)
for i = 1:length(Reb)
   writeVideo(v, MOV(i)) 
end
close(v)