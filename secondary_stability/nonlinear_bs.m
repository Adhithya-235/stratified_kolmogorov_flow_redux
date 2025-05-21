function [Xihat, Bhat, Psihat, Lx, Lz, Nx, Nz] = nonlinear_bs(Reb, rescale)

% Construct nonlinear basic state vector in spectral space. 

%% LOAD DATA FROM FILE

datafile = sprintf('import/ECS_real_field_for_Reb=%d.mat', Reb);
nlbs     = load(datafile);

%% GET ORIGINAL Nx, Nz

Nx = length(nlbs.x);
Nz = length(nlbs.z);

%% PHYSICAL SPACE CONSTRUCTION

Xi  = interpft(interpft(nlbs.omega, Nz*rescale, 1), Nx*rescale, 2);
Psi = interpft(interpft(nlbs.psi, Nz*rescale, 1), Nx*rescale, 2);
B   = interpft(interpft(nlbs.Buoyancy, Nz*rescale, 1), Nx*rescale, 2); 

%% FOURIER SPACE

Xihat  = fft2(Xi);
Psihat = fft2(Psi);
Bhat   = fft2(B);

%% UNWRAP

Xihat  = Xihat(:);
Psihat = Psihat(:);
Bhat   = Bhat(:);

%% GET Lx AND Lz

Lx = 2*nlbs.x(end) - nlbs.x(end-1);
Lz = 2*nlbs.z(end) - nlbs.z(end-1);

%% GET RESCALED Nx, Nz

Nx = Nx*rescale;
Nz = Nz*rescale;

end