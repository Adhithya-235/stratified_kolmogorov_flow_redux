function check_fftSymmetry(F, Nx, Nz)

%% BLOCKWISE DECOMPOSITION

kz01 = F(1, 2:0.5*Nx); kz02 = F(1, (0.5*Nx+2):end);  
kx01 = F(2:0.5*Nz, 1); kx02 = F((0.5*Nz+2):end, 1);
A    = F(2:0.5*Nz, 2:0.5*Nx);
B    = F(2:0.5*Nz, (0.5*Nx+2):end);
C    = F((0.5*Nz+2):end, 2:0.5*Nx);
D    = F((0.5*Nz+2):end, (0.5*Nx+2):end);

%% CHECKS

checkkz0 = max(abs(kz02 - fliplr(conj(kz01))));
checkkx0 = max(abs(kx02 - flipud(conj(kx01))));
checkDA  = max(max(abs(D - rot90(conj(A),2))));
checkBC  = max(max(abs(C - rot90(conj(B),2))));

% disp([checkkz0, checkkx0, checkDA, checkBC])

%% RESULT

cond = abs(checkkz0) + abs(checkkx0) + abs(checkDA) + abs(checkBC) == 0;
msg  = 'Hermitian symmetry violated.';
assert(cond, msg)

end