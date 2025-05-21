function J = jacobn(qb, q, Kx, Kz, alpha, beta, Lx, Lz)
% Applies the jacobian operator to the unwrapped 2D fields qbhat, qhat 
% and returns the unwrapped field Jhat. All fields are given in 
% spectral space. qbhat is a basic state field, alpha and beta 
% are Floquet modifiers. The jacobian is given by
% 
%   J(Q, q) = dz(Q) dx(q) - dx(Q) dz(q)

% disp('Jac Call')
%% GET NO. OF GRID POINTS IN x AND z

Nx = size(Kx,2);
Nz = size(Kx,1);

%% CALCULATE DERIVATIVES IN SPECTRAL SPACE

dx_qb = diffn(qb, Kx);
dz_qb = diffn(qb, Kz);
dx_q  = diffn(q, Kx, Lx, alpha);
dz_q  = diffn(q, Kz, Lz, beta);

%% RESHAPE DERIVATIVES

dx_Qb = reshape(dx_qb, Nz, Nx);
dz_Qb = reshape(dz_qb, Nz, Nx);
dx_Q  = reshape(dx_q, Nz, Nx);
dz_Q  = reshape(dz_q, Nz, Nx);

%% COMPUTE JACOBIAN IN GRID SPACE

J = (dz_Qb.*dx_Q) - (dx_Qb.*dz_Q); 

%% UNWRAP

J = J(:);

end