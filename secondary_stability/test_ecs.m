clear
close all
clc

%% DIARY

diary('logfile001.txt')

%% PHYSICAL PARAMETERS

Reb   = 12;                  % Buoyancy Reynolds
Pr    = 1;                  % Prandtl
Fr    = 0.01;               % Froude
alpha = 0.5;                % x Floquet Modifier
beta  = 0;                  % z Floquet Modifier

%% EIG PARAMETERS

numEig      = 30;
residualtol = 1e-6;
maxit       = 2000;
outputflag  = 1;

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
