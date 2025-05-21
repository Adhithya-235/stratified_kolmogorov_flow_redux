function Aq = boussinesq(q, Xi, B, Psi, Kx, Kz, alpha, beta,...
    Lx, Lz, Fr, Pr, Reb)
% Computes the action of the linear Boussinesq operator upon the state
% vector qhat = (xihat, bhat)'. The inputs are
% Xihat, Bhat, Psihat ---> Basic state vorticity, buoyancy, streamfunction coeffs.
% Kx, Kz ---> Wavenumber matrices/meshes.
% alpha, beta ---> x and z Floquet modifiers. 
% Lx, Lz ---> x and z base state periods.
% Fr, Pr, Reb ---> Froude, Prandtl, Buoyancy Reynolds Numbers. 

% disp('Boussinesq start')

%% DECOMPOSE STATE VECTOR

[xi, b] = decompose_eigenfunction(q);

%% COMPUTE STREAMFUNCTION

psi = inverse_laplacian(xi, Kx, Kz, Fr, alpha, beta, Lx, Lz); 

%% VORTICITY EQUATION TERMS

Dxi      = laplacian(xi, Kx, Kz, Fr, alpha, beta, Lx, Lz);
dxb      = diffn(b, Kx, Lx, alpha);
J_Psi_xi = jacobn(Psi, xi, Kx, Kz, alpha, beta, Lx, Lz);
J_Xi_psi = jacobn(Xi, psi, Kx, Kz, alpha, beta, Lx, Lz);

A1_xi    = (1/Reb)*Dxi - dxb - J_Psi_xi + J_Xi_psi;

%% BUOUYANCY EQUATION TERMS

Db       = laplacian(b, Kx, Kz, Fr, alpha, beta, Lx, Lz);
dxpsi    = diffn(psi, Kx, Lx, alpha);
J_Psi_b  = jacobn(Psi, b, Kx, Kz, alpha, beta, Lx, Lz);
J_B_psi  = jacobn(B, psi, Kx, Kz, alpha, beta, Lx, Lz);

A2_b     = (1/(Pr*Reb))*Db + dxpsi - J_Psi_b + J_B_psi;

%% ASSEMBLE RESULT

Aq = [A1_xi; A2_b];

end