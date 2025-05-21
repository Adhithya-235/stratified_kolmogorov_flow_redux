function Dqhat = laplacian(qhat, Kx, Kz, Fr, alpha, beta, Lx, Lz)
% Applies the Laplacian operator to the unwrapped 2D field qhat 
% and returns the unwrapped field Dqhat. The field qhat is given in 
% spectral space, and so is Dqhat. Fr is the Froude Number; alpha and beta 
% are Floquet modifiers. 

%% FUNDAMENTAL WAVENUMBERS

kx = 2*pi/Lx;
kz = 2*pi/Lz;

%% GET NO. OF GRID POINTS IN x AND z

Nx = size(Kx,2);
Nz = size(Kx,1);

%% FORM IDENTITY MATRIX

I = ones(Nz,Nx);

%% COMPUTE LAPLACIAN

Qhat  = reshape(qhat, Nz, Nx); 
DQhat = -((Fr.^2).*((Kx + kx*alpha*I).^2) + (Kz + kz*beta*I).^2).*Qhat;
Dqhat = DQhat(:);

end