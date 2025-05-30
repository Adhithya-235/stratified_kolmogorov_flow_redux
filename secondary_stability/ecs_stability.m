% ecs_stability.m
%
% This script performs secondary stability analysis on a given ECS for
% specified parameters. It expects the following variables to be defined
% in the MATLAB workspace prior to execution:
%
%   Reb         - Buoyancy Reynolds number
%   Pr          - Prandtl number
%   Fr          - Froude number
%   alpha       - Floquet exponent in x
%   beta        - Floquet exponent in z
%   numEig      - Number of eigenvalues to compute
%   residualtol - Arnoldi residual tolerance
%   maxit       - Maximum Arnoldi iterations
%   outputflag  - Verbosity flag (1 = verbose)
%
% Note: This script assumes that the commands `clear`, `close all`, and `clc`
% have already been issued in the MATLAB command sent via `-batch`.
%
% This script is designed to be called from the shell using:
%   matlab -batch "clear; close all; clc; Reb=12; Pr=1; Fr=0.01; alpha=0.5; beta=0; numEig=30; residualtol=1e-6; maxit=2000; outputflag=1; ecs_stability"

%% CHECK FOR REQUIRED INPUTS

required_vars = {'Reb','Pr','Fr','alpha','beta','numEig','residualtol','maxit','outputflag'};
for k = 1:length(required_vars)
    if ~exist(required_vars{k}, 'var')
        error('Missing required variable: %s', required_vars{k});
    end
end

%% SAVE FOLDER

solDir = fullfile("solutions", sprintf("Reb%.2f_alpha%.2f", Reb, alpha));
if ~exist(solDir, 'dir')
    mkdir(solDir);
end

%% LOGGING ON

logFile = fullfile(solDir, sprintf('ecs_convergence_log_Reb%.2f_alpha%.2f.txt', Reb, alpha));
if exist(logFile, 'file')
    delete(logFile);  % delete old log file
end
diary(logFile);
diary on;

%% MAIN SECTION START

try

    [eigvals, eigvecs, dom_mode, meta] = compute_converged_spectrum(Reb, Pr,...
        Fr, alpha, beta, numEig, residualtol, maxit);

    %% LIST EIGENVALUES

    disp('Computed eigenvalues:');
    disp(eigvals);

    %% PLOT DOMINANT EIGENFUNCTION
    
    fprintf("Plotting dominant eigenvectors...\n");
    [Xp, Zp] = meshgrid(dom_mode.xp, dom_mode.zp);
    dommode = combinedEigenfunctionPlot(Xp/Fr, Zp, dom_mode.Xi, dom_mode.Psi, dom_mode.B);

    %% PRINT TIMES

    fprintf('Completed successfully. Iteration time: %.4f seconds\n', meta.itertime);

    %% SAVE SOLUTION TO .MAT FILE

    solFile = fullfile(solDir, sprintf("spectrum_Reb%.2f_alpha%.2f.mat", Reb, alpha));
    save(solFile, "eigvals", "eigvecs", "dom_mode", "meta");

    %% PLOT EIGENVALUES

    fprintf("Plotting eigenvalue spectrum...\n");

    espec = figure('Color', 'w', 'WindowState','maximized');
    theta = linspace(0, 2*pi, 500);
    plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.5); hold on; % unit circle
    scatter(real(eigvals), imag(eigvals), 40, 'filled');
    xlabel('$\mathrm{Fr} \cdot \mathrm{Re}(\sigma)$', 'Interpreter', 'latex', 'FontSize', 30);
    ylabel('$\mathrm{Fr} \cdot \mathrm{Im}(\sigma)$', 'Interpreter', 'latex', 'FontSize', 30);
    xlim([-2, 2])
    ylim([-2, 2])
    axis square;
    grid on;
    box on
    set(gca, 'FontSize', 30, 'LineWidth', 3);

    %% SAVE EIGENVALUE PLOT

    figPath = fullfile(solDir, sprintf("spectrum_Reb%.2f_alpha%.2f.fig", Reb, alpha));
    pngPath = fullfile(solDir, sprintf("spectrum_Reb%.2f_alpha%.2f.png", Reb, alpha));
    savefig(espec, figPath);
    exportgraphics(espec, pngPath, 'Resolution', 500);

    %% SAVE DOMINANT MODE PLOT

    fprintf("Saving dominant mode plots...\n");
    domFigPath = fullfile(solDir, sprintf("dommode_Reb%.2f_alpha%.2f.fig", Reb, alpha));
    domPngPath = fullfile(solDir, sprintf("dommode_Reb%.2f_alpha%.2f.png", Reb, alpha));
    savefig(dommode, domFigPath);
    exportgraphics(dommode, domPngPath, 'Resolution', 500);

%% CATCH EXCEPTION

catch ME
    
    diary off;
    rethrow(ME);

end

%% NORMAL END

diary off
