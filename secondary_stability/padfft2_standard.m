function Fpad = padfft2_standard(F, Nz_new, Nx_new)
    % Zero-pad a 2D FFT array (standard ordering) from size [Nz, Nx] to [Nz_new, Nx_new]

    [Nz, Nx] = size(F);
    Fpad = zeros(Nz_new, Nx_new);

    % Define FFT frequency indices (standard order)
    kx_c = [0:Nx/2-1, -Nx/2:-1];
    kz_c = [0:Nz/2-1, -Nz/2:-1];
    kx_f = [0:Nx_new/2-1, -Nx_new/2:-1];
    kz_f = [0:Nz_new/2-1, -Nz_new/2:-1];

    % Find where each coarse index sits in the fine grid
    [~, kx_idx] = ismember(kx_c, kx_f);
    [~, kz_idx] = ismember(kz_c, kz_f);

    % Insert coarse values into fine grid
    for iz = 1:Nz
        for ix = 1:Nx
            i_fine = kz_idx(iz);
            j_fine = kx_idx(ix);
            if i_fine > 0 && j_fine > 0
                Fpad(i_fine, j_fine) = F(iz, ix);
            end
        end
    end
end