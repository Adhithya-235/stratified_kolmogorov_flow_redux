function Aqhat = boussinesq(qhat, Xihat, Bhat, Psihat, Kx, Kz, alpha, beta,...
    Lx, Lz, Fr, Pr, Reb)
% Computes the action of the linear Boussinesq operator upon the state
% vector qhat = (xihat, bhat)'. The inputs are
% Xihat, Bhat, Psihat ---> Basic state vorticity, buoyancy, streamfunction coeffs.
% Kx, Kz ---> Wavenumber matrices/meshes.
% alpha, beta ---> x and z Floquet modifiers. 
% Lx, Lz ---> x and z base state periods.
% Fr, Pr, Reb ---> Froude, Prandtl, Buoyancy Reynolds Numbers. 

%% DECOMPOSE STATE VECTOR

[xihat, bhat] = decompose_eigenfunction(qhat);

%% COMPUTE STREAMFUNCTION

psihat = inverse_laplacian(xihat, Kx, Kz, Fr, alpha, beta, Lx, Lz); 

%% VORTICITY EQUATION TERMS

Dxihat      = laplacian(xihat, Kx, Kz, Fr, alpha, beta, Lx, Lz);
dxbhat      = diffn(bhat, Kx, Lx, alpha);
Jhat_Psi_xi = jacobn(Psihat, xihat, Kx, Kz, alpha, beta, Lx, Lz);
Jhat_Xi_psi = jacobn(Xihat, psihat, Kx, Kz, alpha, beta, Lx, Lz);

A1_xihat    = (1/Reb)*Dxihat - dxbhat - Jhat_Psi_xi + Jhat_Xi_psi;

%% BUOUYANCY EQUATION TERMS

Dbhat       = laplacian(bhat, Kx, Kz, Fr, alpha, beta, Lx, Lz);
dxpsihat    = diffn(psihat, Kx, Lx, alpha);
Jhat_Psi_b  = jacobn(Psihat, bhat, Kx, Kz, alpha, beta, Lx, Lz);
Jhat_B_psi  = jacobn(Bhat, psihat, Kx, Kz, alpha, beta, Lx, Lz);

A2_bhat     = (1/(Pr*Reb))*Dbhat + dxpsihat - Jhat_Psi_b + Jhat_B_psi;

%% ASSEMBLE RESULT

Aqhat = [A1_xihat; A2_bhat];

end