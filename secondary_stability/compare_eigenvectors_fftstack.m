function vecErr = compare_eigenvectors_fftstack(v_coarse, v_fine, Nx_coarse, Nz_coarse, Nx_fine, Nz_fine)
    % Split and reshape coarse
    MN_coarse = Nx_coarse*Nz_coarse;
    xi_coarse = reshape(v_coarse(1:MN_coarse), Nz_coarse, Nx_coarse);
    b_coarse  = reshape(v_coarse(MN_coarse+1:2*MN_coarse), Nz_coarse, Nx_coarse);

    % Split and reshape fine
    MN_fine = Nx_fine*Nz_fine;
    xi_fine = reshape(v_fine(1:MN_fine), Nz_fine, Nx_fine);
    b_fine  = reshape(v_fine(MN_fine+1:2*MN_fine), Nz_fine, Nx_fine);

    % Pad coarse to fine size
    xi_interp = padfft2_standard(xi_coarse, Nz_fine, Nx_fine);
    b_interp  = padfft2_standard(b_coarse, Nz_fine, Nx_fine);

    % Stack and normalize
    v1 = [xi_interp(:); b_interp(:)];
    v2 = [xi_fine(:);   b_fine(:)];

    v1 = v1 / norm(v1);
    v2 = v2 / norm(v2);

    % Compute 1 - cosine similarity
    vecErr = 1 - abs(v1' * v2);
end