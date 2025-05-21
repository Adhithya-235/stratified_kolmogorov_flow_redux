function [Xi, B, Psi] = laminar_bs(Z, m)

% Construct laminar basic state vector in grid space. 

%% PHYSICAL SPACE CONSTRUCTION

Xi  = m*cos(m*Z);
Psi = -(1/m)*cos(m*Z);
B   = zeros(size(Z)); 

%% UNWRAP

Xi  = Xi(:);
Psi = Psi(:);
B   = B(:);

end