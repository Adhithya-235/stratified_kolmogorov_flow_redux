clear
close all
clc

%% ADD UTILITIES TO PATH

addpath('../../utility_belt');
addpath('..'); 

%% PLOT AESTHETICS

fs   = 24;
lw   = 3;

%% PARAMETERS

Reb     = [10, 11, 12, 13, 15, 20, 30, 50, 55, 60, 85, 90, 100];
Fr      = 0.01;
alpha   = 0;
beta    = 0;

%% GET FILE NAMES

SolnFilePath = cell(length(Reb),1);
for fileidx = 1:length(Reb)
    ParentDir             = sprintf('../secondary_stability_solutions_maxit10000');
    subDir                = sprintf('Reb%.2f_alpha%.2f',Reb(fileidx), alpha);
    SolnFile              = sprintf('spectrum_Reb%.2f_alpha%.2f.mat',Reb(fileidx), alpha);
    SolnFilePath{fileidx} = sprintf('%s/%s/%s',ParentDir,subDir,SolnFile);
end

%% PREALLOCATE SOLUTION

nReb         = length(Reb);
tracked_eigs = NaN(nReb, 1);
tracked_idx  = NaN(nReb, 1);

%% GET TARGET MODE

target          = load(SolnFilePath{1});
ref_idx         = 1;
eigval_prev     = target.eigvals(ref_idx);
eigvec_prev     = target.eigvecs(:, ref_idx);

tracked_eigs(1) = eigval_prev;
tracked_idx(1)  = ref_idx;
Nx_coarse       = target.meta.resolution(1);
Nz_coarse       = target.meta.resolution(2);

%% START LOGGING 

fprintf('Tracking mode at Reb = %.2f with sigma = %.4f.\n', Reb(1), eigval_prev)

%% LOOP THROUGH THE DATA

for idx = 2:nReb

    %% LOAD FILE AND GET EIGENVALS, EIGENVECS, AND META

    candidates = load(SolnFilePath{idx});
    sigma      = candidates.eigvals;
    V          = candidates.eigvecs;
    Nx_fine    = candidates.meta.resolution(1);
    Nz_fine    = candidates.meta.resolution(2);

    %% PAD REFERENCE EIGENVECTOR IF DIMENSIONS DON'T MATCH 

    if (Nx_coarse ~= Nx_fine) || (Nz_coarse ~= Nz_fine)
        v_ref = pad_coarse_evec(eigvec_prev, Nx_coarse, Nz_coarse, Nx_fine, Nz_fine);
    else
        v_ref = eigvec_prev;
    end

    %% COMPUTE OVERLAP

    overlaps = abs(V'*v_ref);

    %% CHOOSE BEST MATCH

    [best_overlap, best_match_idx] = max(overlaps);
    eigval_prev                    = sigma(best_match_idx);
    eigvec_prev                    = V(:, best_match_idx);
    Nx_coarse                      = Nx_fine;
    Nz_coarse                      = Nz_fine;

    %% SAVE IN SOLUTION ARRAY

    tracked_eigs(idx) = eigval_prev;
    tracked_idx(idx)  = best_match_idx;

    %% LOG RESULT

    fprintf('Match found at Reb = %.2f with sigma = %.4f and overlap %.3f.\n', Reb(idx), eigval_prev, best_overlap)


end

%% PLOT TRACKED EIGENVALUES

trackedevals = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(Reb, real(tracked_eigs), '-o', 'linewidth', lw)
xlabel('$\mathrm{Re}_b$', 'interpreter', 'latex')
ylabel('$\mathrm{Fr}\,\mathcal{R}\left[\sigma\right]$', 'interpreter', 'latex')
set(gca, 'fontsize', fs)
xlim([1, Reb(end)])
axis square
grid on
box on
set(gca, 'linewidth', lw, 'fontsize', fs, 'XScale', 'linear', 'YScale', 'linear', 'ticklabelinterpreter', 'latex', 'GridLineStyle', ':')
drawnow
plotname = sprintf('StationaryInstabGrowth.pdf');
exportgraphics(trackedevals, plotname, 'ContentType', 'vector', 'Resolution', 500);

%% EIGENVECTOR PADDING FUNCTION

function v_interp = pad_coarse_evec(v_coarse, Nx_coarse, Nz_coarse, Nx_fine, Nz_fine)
    
    % Split and reshape coarse
    MN_coarse = Nx_coarse*Nz_coarse;
    xi_coarse = reshape(v_coarse(1:MN_coarse), Nz_coarse, Nx_coarse);
    b_coarse  = reshape(v_coarse(MN_coarse+1:2*MN_coarse), Nz_coarse, Nx_coarse);

    % Pad coarse to fine size
    xi_interp = padfft2_standard(xi_coarse, Nz_fine, Nx_fine);
    b_interp  = padfft2_standard(b_coarse, Nz_fine, Nx_fine);

    % Stack and normalize
    v_interp = [xi_interp(:); b_interp(:)];
    v_interp = v_interp / norm(v_interp);

end