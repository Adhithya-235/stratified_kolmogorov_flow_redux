function FH = enforce_fftSymmetry(F, Nx, Nz)

%% BLOCKWISE INDEX DECOMPOSITION

kz01 = {1, 2:0.5*Nx}; kz02 = {1, (0.5*Nx+2):Nx};               % kz = 0
kzN1 = {0.5*Nz+1, 2:0.5*Nx}; kzN2 = {0.5*Nz+1, (0.5*Nx+2):Nx}; % kz = Nz/2
kx01 = {2:0.5*Nz, 1}; kx02 = {(0.5*Nz+2):Nz, 1};               % kx = 0  
kxN1 = {2:0.5*Nz, 0.5*Nx+1}; kxN2 = {(0.5*Nz+2):Nz, 0.5*Nx+1}; % kx = Nx/2
A    = {2:0.5*Nz, 2:0.5*Nx};                                   % UL quadrant
B    = {2:0.5*Nz, (0.5*Nx+2):Nx};                              % UR quadrant
C    = {(0.5*Nz+2):Nz, 2:0.5*Nx};                              % LL quadrant
D    = {(0.5*Nz+2):Nz, (0.5*Nx+2):Nz};                         % LR quadrant

%% APPLY CONJUGATE RULES TO ZERO AND NYQUIST MODES

F(kz02{:}) = fliplr(conj(F(kz01{:}))); 
F(kzN2{:}) = fliplr(conj(F(kzN1{:})));
F(kx02{:}) = flipud(conj(F(kx01{:})));
F(kxN2{:}) = flipud(conj(F(kxN1{:})));

%% APPLY CONJUGATE RULES TO REMAINING QUADRANTS

F(D{:}) = rot90(conj(F(A{:})), 2);
F(C{:}) = rot90(conj(F(B{:})), 2);

%% SAVE IN DIFFERENT MATRIX

FH = F;

end