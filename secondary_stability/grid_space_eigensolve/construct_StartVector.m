function q0hat = construct_StartVector(Nx, Nz); 
% Construct matrices that would be produced by ifft-ing a real field to use
% to start the Arnoldi routine. 

%% GET RANDOM FIELDS

Xi0 = randn(Nz,Nx);
B0  = randn(Nz,Nx); 

%% CHECK HERMITIAN 

check_fftSymmetry(fft2(Xi0), Nx, Nz);
check_fftSymmetry(fft2(B0), Nx, Nz);

%% COMPOSE

q0hat = [Xi0(:); B0(:)];

end