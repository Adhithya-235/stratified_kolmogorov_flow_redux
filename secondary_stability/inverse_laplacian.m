function qhat = inverse_laplacian(Dqhat, Kx, Kz, Fr, alpha, beta, Lx, Lz)
% Applies the inverse Laplacian operator to the unwrapped 2D field Dqhat 
% and returns the unwrapped field qhat. The field Dqhat is given in 
% spectral space, and so is qhat. Fr is the Froude Number; alpha and beta 
% are Floquet modifiers. 

%% FUNDAMENTAL WAVENUMBERS

kx = 2*pi/Lx;
kz = 2*pi/Lz;

%% GET NO. OF GRID POINTS IN x AND z

Nx = size(Kx,2);
Nz = size(Kx,1);

%% FORM IDENTITY MATRIX

I = ones(Nz,Nx);

%% DEFINE LAPLACIAN 

L0 = -(Kx.^2 + Kz.^2);
L  = -((Fr.^2).*((Kx + kx*alpha*I).^2) + (Kz + kz*beta*I).^2);

%% COMPUTE INVERSE LAPLACIAN

DQhat       = reshape(Dqhat, Nz, Nx); 
Qhat        = zeros(size(DQhat));
Qhat(L0~=0) = DQhat(L0~=0)./L(L0~=0);
Qhat(L0==0) = 0;
qhat = Qhat(:);

end