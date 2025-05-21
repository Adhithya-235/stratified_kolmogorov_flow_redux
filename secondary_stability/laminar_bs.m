function [Xihat, Bhat, Psihat] = laminar_bs(Z, m)

% Construct laminar basic state vector in spectral space. 

%% PHYSICAL SPACE CONSTRUCTION

Xi  = m*cos(m*Z);
Psi = -(1/m)*cos(m*Z);
B   = zeros(size(Z)); 

%% FOURIER SPACE

Xihat  = fft2(Xi);
Psihat = fft2(Psi);
Bhat   = fft2(B);

%% UNWRAP

Xihat  = Xihat(:);
Psihat = Psihat(:);
Bhat   = Bhat(:);

end