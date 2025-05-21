function q0hat = construct_StartVector(Nx, Nz); 
% Construct matrices that would be produced by ifft-ing a real field to use
% to start the Arnoldi routine. 

%% GET RANDOM FIELDS

Xi0 = randn(Nz,Nx);
B0  = randn(Nz,Nx); 

%% DO AN FFT

Xi0hat = fft2(Xi0);
B0hat  = fft2(B0); 

%% UNWRAP

xi0hat = Xi0hat(:);
b0hat  = B0hat(:);

%% COMPOSE

q0hat = [xi0hat; b0hat];

end