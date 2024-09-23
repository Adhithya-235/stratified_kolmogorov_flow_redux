function [Ub, Ub2, Bb1, z] = construct_base_state(N, H, m)

% Construct a sinusoidal base state for U with vertical wavenumber m on the
% chosen grid with N, H.

%% GET GRID AND MATRICES TO CONSTRUCT BASIC STATE

[x, D2] = difmat(N,2); 
scale = H/(2*pi);
z = scale*x;
D2 = (1/(scale^2))*D2;

%% CONSTRUCT BASIC STATE

Ub  = sin(m*z);
Ub2 = D2*Ub;
Bb1 = 0*Ub;

end