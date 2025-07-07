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

Reb            = [1:20, 30:10:100];
Fr             = 0.01;
alpha          = 0;
beta           = 0;
Reb_start      = 18;
ref_idx        = 1;
overlap_thresh = 0.5;
output_file    = 'oscillatory_mode_1_data.mat';

%% GET FILE NAMES

SolnFilePath = cell(length(Reb),1);
for fileidx = 1:length(Reb)
    ParentDir             = sprintf('../solutions_branch_1_Pr_1');
    subDir                = sprintf('Reb%.2f_alpha%.2f',Reb(fileidx), alpha);
    SolnFile              = sprintf('spectrum_Reb%.2f_alpha%.2f.mat',Reb(fileidx), alpha);
    SolnFilePath{fileidx} = sprintf('%s/%s/%s',ParentDir,subDir,SolnFile);
end

%% PREALLOCATE SOLUTION

nReb         = length(Reb);
tracked_eigs = NaN(nReb, 1);
tracked_idx  = NaN(nReb, 1);
tracked_vecs = cell(nReb, 1);

%% IDENTIFY START INDEX

start_idx = find(Reb == Reb_start);
if isempty(start_idx)
    error('Reb_start = %.2f not found in Reb array.', Reb_start)
end

%% GET TARGET MODE

target          = load(SolnFilePath{start_idx});
eigval_prev     = target.eigvals(ref_idx);
eigvec_prev     = target.eigvecs(:, ref_idx);
Nx_coarse       = target.meta.resolution(1);
Nz_coarse       = target.meta.resolution(2);

tracked_eigs(start_idx) = eigval_prev;
tracked_idx(start_idx)  = ref_idx;
tracked_vecs{start_idx} = eigvec_prev;

%% START LOGGING 

fprintf('Tracking mode at Reb = %.2f with sigma = %.4f.\n', Reb(start_idx), eigval_prev)

%% FORWARD TRACKING

if start_idx < nReb
    
    % re-initialize the reference grid size at the start point
    data0        = load(SolnFilePath{start_idx});
    Nx_coarse    = data0.meta.resolution(1);
    Nz_coarse    = data0.meta.resolution(2);
    eigvec_prev  = tracked_vecs{start_idx};

    for idx = start_idx+1:nReb
    
        [eigval_prev, eigvec_prev, Nx_coarse, Nz_coarse, best_idx, overlap] = ...
            track_one_step(SolnFilePath{idx}, eigvec_prev, Nx_coarse, Nz_coarse);

        if overlap < overlap_thresh
            fprintf('Reb=%.2f: overlap %.3f < %.2f → discarded\n', Reb(idx), overlap, overlap_thresh);
            continue
        end

        tracked_eigs(idx)  = eigval_prev;
        tracked_idx(idx)   = best_idx;
        tracked_vecs{idx}  = eigvec_prev;

        fprintf('Match at Reb = %.2f | σ = %.4f | overlap = %.3f\n', ...
            Reb(idx), eigval_prev, overlap);
    
    end

end

%% BACKWARD TRACKING

if start_idx > 1

    % again re-init grid size at the start point
    data0        = load(SolnFilePath{start_idx});
    Nx_coarse    = data0.meta.resolution(1);
    Nz_coarse    = data0.meta.resolution(2);
    eigvec_prev  = tracked_vecs{start_idx};

    for idx = start_idx-1:-1:1
        
        [eigval_prev, eigvec_prev, Nx_coarse, Nz_coarse, best_idx, overlap] = ...
            track_one_step(SolnFilePath{idx}, eigvec_prev, Nx_coarse, Nz_coarse);

        if overlap < overlap_thresh
            fprintf('Reb=%.2f: overlap %.3f < %.2f → discarded\n', Reb(idx), overlap, overlap_thresh);
            continue
        end

        tracked_eigs(idx)  = eigval_prev;
        tracked_idx(idx)   = best_idx;
        tracked_vecs{idx}  = eigvec_prev;

        fprintf('Match at Reb = %.2f | σ = %.4f | overlap = %.3f\n', ...
            Reb(idx), eigval_prev, overlap);

    end

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
set(gca, 'linewidth', lw, 'fontsize', fs, 'XScale', 'log', 'YScale', 'linear', 'ticklabelinterpreter', 'latex', 'GridLineStyle', ':')
drawnow
plotname = sprintf('Oscillatory1InstabGrowth.png');
exportgraphics(trackedevals, plotname, 'ContentType', 'vector', 'Resolution', 500);

%% SAVE TRACKED DATA

save(output_file, 'Reb', 'tracked_eigs', 'tracked_idx', 'tracked_vecs');
fprintf('Saved tracked eigenvalues and vectors to %s.\n', output_file);


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

%% TRACK ONE STEP FUNCTION

function [eigval_new, eigvec_new, Nx_fine, Nz_fine, best_idx, overlap] = ...
    track_one_step(file_path, v_ref, Nx_coarse, Nz_coarse)

    data     = load(file_path);
    sigma    = data.eigvals;
    V        = data.eigvecs;
    Nx_fine  = data.meta.resolution(1);
    Nz_fine  = data.meta.resolution(2);

    if (Nx_coarse ~= Nx_fine) || (Nz_coarse ~= Nz_fine)
        v_ref = pad_coarse_evec(v_ref, Nx_coarse, Nz_coarse, Nx_fine, Nz_fine);
    end

    overlaps     = abs(V' * v_ref);
    invalid = imag(sigma) < -1e-5;
    overlaps(invalid) = -Inf;
    [overlap, best_idx] = max(overlaps);
    eigval_new   = sigma(best_idx);
    eigvec_new   = V(:, best_idx);

end
