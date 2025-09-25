% jacobn_equivalence_test_128.m
clear 
close all
clc

rng(1);

% --- grid & wavenumbers (unshifted FFT layout) ---
Nx=128; Nz=128; Lx=2*pi; Lz=2*pi; beta=0;
mx = 0:Nx-1; mx(mx> Nx/2)=mx(mx> Nx/2)-Nx; kx=(2*pi/Lx)*mx;
mz = 0:Nz-1; mz(mz> Nz/2)=mz(mz> Nz/2)-Nz; kz=(2*pi/Lz)*mz;
[KX,KZ] = meshgrid(kx,kz);     % Nz x Nx

% --- random test spectra on the Lx-tile (base Q, envelope q) ---
Qhat = randn(Nz,Nx)+1i*randn(Nz,Nx);
qhat = randn(Nz,Nx)+1i*randn(Nz,Nx);

% --- helper: 3/2 padding/truncation (centered) ---
Nxp=floor(3*Nx/2); Nzp=floor(3*Nz/2); sx=floor((Nxp-Nx)/2); sy=floor((Nzp-Nz)/2);

% (inline versions using sy,sx)
pad3 = @(H) ifftshift(padarray(fftshift(H),[sy sx],0,'both'));
% trunc3 = @(Hp) ifftshift(fftshift(Hp)(1+sy:sy+Nz, 1+sx:sx+Nx)); %#ok<NASGU> % (kept for clarity)

% --- padded wavenumbers ---
mxp=0:Nxp-1; mxp(mxp> Nxp/2)=mxp(mxp> Nxp/2)-Nxp; kxp=(2*pi/Lx)*mxp;
mzp=0:Nzp-1; mzp(mzp> Nzp/2)=mzp(mzp> Nzp/2)-Nzp; kzp=(2*pi/Lz)*mzp;
[KXP,KZP] = meshgrid(kxp,kzp);   % Nzp x Nxp
xp = (0:Nx-1)*(Lx/Nx);       % padded physical x-grid
xpp = (0:Nxp-1)*(Lx/Nxp);  
phase_up = exp(1i*(pi/Lx)*xp); % for alpha=1/2; will set alpha below

for alpha = [0, 0.5]
  % ---------------------- Path A: your jacobn (on Lx tile) ----------------------
  J1 = jacobn(Qhat(:), qhat(:), KX, KZ, alpha, 0, Lx, Lz);
  J1 = reshape(J1, Nz, Nx);

  % ---------------- Path B: Bloch → plain J → demodulate (dealiased) -----------
  % pad spectra to padded grid
  Qhat_p = Qhat;
  qhat_p = qhat;

  % envelope on padded grid, then Bloch-modulate in physical space
  q_env_up = ifft2(qhat_p);
  phase = exp(1i*(alpha*2*pi/Lx)*xp); phase = repmat(phase, Nz, 1);
  phasep = exp(1i*(alpha*2*pi/Lx)*xpp); phasep = repmat(phasep, Nzp, 1);
  q_phys_up = q_env_up .* phase;

  % base derivatives on padded grid (plain)
  Qx_up = ifft2(pad3(1i*KX .* Qhat_p));
  Qz_up = ifft2(pad3(1i*KZ .* Qhat_p));

  % perturbation derivatives of the physical (Bloch) field (plain)
  qphys_hat_up = fft2(q_phys_up);
  qx_up = ifft2(pad3(1i*KX.* qphys_hat_up));
  qz_up = ifft2(pad3(1i*KZ.* qphys_hat_up));

  % physical-space Jacobian (padded), then demodulate back to envelope
  J_phys_up = Qz_up .* qx_up - Qx_up .* qz_up;
  J_env_up  = J_phys_up ./ phasep;

  % FFT and truncate back to Lx×Lz
  J2_up_hat = fft2(J_env_up);
  C = fftshift(J2_up_hat);
  Csmall = C(1+sy:sy+Nz, 1+sx:sx+Nx);
  J2 = ifftshift(Csmall);     % Nz x Nx (unshifted), envelope frame

  % ------------------------------ compare --------------------------------------
  rel_err = norm(J1 - J2, 'fro') / max(norm(J2,'fro'), eps);
  fprintf('alpha = %.1f :  ||J_Floquet - J_Bloch|| / ||J_Bloch|| = %.3e\n', alpha, rel_err);
end
