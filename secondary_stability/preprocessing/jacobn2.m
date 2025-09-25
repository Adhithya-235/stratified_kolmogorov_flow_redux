function Jhat = jacobn2(qbhat, qhat, Kx, Kz, alpha, beta, Lx, Lz)
% Applies the jacobian operator to the unwrapped 2D fields qbhat, qhat 
% and returns the unwrapped field Jhat. All fields are given in 
% spectral space. qbhat is a basic state field, alpha and beta 
% are Floquet modifiers. The jacobian is given by
% 
%   J(Q, q) = dz(Q) dx(q) - dx(Q) dz(q)

%% GET NO. OF GRID POINTS IN x AND z

Nx = size(Kx,2);
Nz = size(Kx,1);

%% CALCULATE DERIVATIVES IN SPECTRAL SPACE

dx_qbhat = diffn(qbhat, Kx);
dz_qbhat = diffn(qbhat, Kz);
dx_qhat  = diffn(qhat, Kx, Lx, alpha);
dz_qhat  = diffn(qhat, Kz, Lz, beta);

%% RESHAPE DERIVATIVES

dx_Qbhat = reshape(dx_qbhat, Nz, Nx);
dz_Qbhat = reshape(dz_qbhat, Nz, Nx);
dx_Qhat  = reshape(dx_qhat, Nz, Nx);
dz_Qhat  = reshape(dz_qhat, Nz, Nx);

%% PERFORM IFFTs

dx_Qb = ifft2(dx_Qbhat);
dz_Qb = ifft2(dz_Qbhat);
dx_Q  = ifft2(dx_Qhat);
dz_Q  = ifft2(dz_Qhat);

%% COMPUTE JACOBIAN IN GRID SPACE

J = (dz_Qb.*dx_Q) - (dx_Qb.*dz_Q); 

%% BACK TO FOURIER SPACE AND UNWRAP

Jhat = fft2(J);
Jhat = Jhat(:);

end