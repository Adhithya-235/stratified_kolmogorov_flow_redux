% SCRIPT: plot_selected_eigenfunction_v2.m
%
% This script prompts the user for Reb, alpha, and an eigenfunction index.
% It then loads the corresponding converged eigenpair data, processes the
% selected eigenfunction using process_eigenfunction.m (by passing the
% full V matrix and the selected index), and plots the resulting 
% Xi, B, and Psi fields.
%
% Assumes that 'process_eigenfunction.m' and its dependencies are in the 
% parent directory of this script's location.
%

clear;
clc;
close all; 

%% SCRIPT INITIALIZATION & CONFIGURATION

addpath('../../utility_belt'); % User's added path
addpath('..'); % User's added path
fprintf('Starting script: plot_selected_eigenfunction_v2.m\n');

% Aesthetic configurations
plot_fontsize = 16;
plot_borderwidth = 3;
% Fr scaling assumed for the eigenvalues stored in the converged_...mat files
Fr_assumed_for_eigenvalues = 0.01; 

%% USER INPUTS

Reb_val   = 30;
alpha_val = 0;
eigenfunction_index_to_plot = 3; % Renamed for clarity

if isnan(Reb_val) || isnan(alpha_val) || isnan(eigenfunction_index_to_plot) || ...
   eigenfunction_index_to_plot < 1 || floor(eigenfunction_index_to_plot) ~= eigenfunction_index_to_plot
    error('Invalid input. Reb, alpha must be numbers, and index must be a positive integer.');
end

fprintf('Attempting to plot eigenfunction %d for Reb = %.2f, alpha = %.2f\n', ...
        eigenfunction_index_to_plot, Reb_val, alpha_val);

%% PATH AND FILE NAME CONSTRUCTION
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir) 
    scriptDir = pwd; 
    fprintf('Assuming script is running from: %s (postprocessing folder)\n', scriptDir);
end
solutionsParentDir = fullfile(scriptDir, '..', 'solutions_branch_1_Pr_1');
parentDir = fullfile(scriptDir, '..'); 

Reb_str_fs = sprintf('%.2f', Reb_val);       
alpha_str_fs = sprintf('%.2f', alpha_val);   

dirName = sprintf('Reb%s_alpha%s', Reb_str_fs, alpha_str_fs);
currentSolutionFolderPath = fullfile(solutionsParentDir, dirName);

Reb_str_for_name = strrep(Reb_str_fs, '.', 'p');     
alpha_str_for_name = strrep(alpha_str_fs, '.', 'p'); 
    
convergedMatFilename = sprintf('converged_Reb%s_alpha%s.mat', Reb_str_for_name, alpha_str_for_name);
convergedMatFilePath = fullfile(currentSolutionFolderPath, convergedMatFilename);
dynamicVarName = sprintf('converged_Reb%s_alpha%s', Reb_str_for_name, alpha_str_for_name);

if ~isfile(convergedMatFilePath)
    error('Converged data file not found: %s\nPlease ensure Reb and alpha are correct and isolate_converged_eigenpairs.m has been run.', convergedMatFilePath);
end

%% LOAD CONVERGED DATA
fprintf('Loading converged data from: %s\n', convergedMatFilename);
loadedFile = load(convergedMatFilePath);

if ~isfield(loadedFile, dynamicVarName)
    error('Dynamically named structure ''%s'' not found in MAT file: %s', dynamicVarName, convergedMatFilename);
end
dataStruct = loadedFile.(dynamicVarName);

requiredStructFields = {'converged_eigenvalues', 'converged_eigenvectors', ...
                        'original_resolution', 'original_domainsize', ...
                        'Fr_scaling_of_eigenvalues', 'alpha'};
for f_idx = 1:length(requiredStructFields)
    if ~isfield(dataStruct, requiredStructFields{f_idx})
        error('Loaded data structure from %s is missing field: %s', convergedMatFilename, requiredStructFields{f_idx});
    end
end

if eigenfunction_index_to_plot > size(dataStruct.converged_eigenvectors, 2) || ...
   eigenfunction_index_to_plot > length(dataStruct.converged_eigenvalues)
    error('Selected eigenfunction index (%d) is out of bounds. File contains %d converged eigenpairs.', ...
          eigenfunction_index_to_plot, length(dataStruct.converged_eigenvalues));
end

%% EXTRACT DATA AND PREPARE FOR process_eigenfunction
V_full_converged = dataStruct.converged_eigenvectors; % Full matrix of converged eigenvectors
selected_eigval = dataStruct.converged_eigenvalues(eigenfunction_index_to_plot); % For title

Nx = dataStruct.original_resolution(1);
Nz = dataStruct.original_resolution(2);
Lxp_stored = dataStruct.original_domainsize(1); 
Lzp_stored = dataStruct.original_domainsize(2); 

alpha_for_scaling = dataStruct.alpha; 
scale_factor = alpha_for_scaling;
if abs(scale_factor) == 0
    scale_factor = 1.0;
end
Lx_computational = Lxp_stored * scale_factor;
Lz_computational = Lzp_stored; 

kx_vals = (2*pi/Lx_computational) * [0:(Nx/2-1), (-Nx/2):-1];
kz_vals = (2*pi/Lz_computational) * [0:(Nz/2-1), (-Nz/2):-1];
[Kx_mesh, Kz_mesh] = meshgrid(kx_vals, kz_vals);

Fr_to_pass    = dataStruct.Fr_scaling_of_eigenvalues; 
alpha_to_pass = dataStruct.alpha; 
beta_to_pass  = 0.0; 
Fr            = Fr_to_pass;

fprintf('Processing eigenfunction using parameters:\n');
fprintf('  Nx=%d, Nz=%d, Lx_comp=%.2f, Lz_comp=%.2f\n', Nx, Nz, Lx_computational, Lz_computational);
fprintf('  Fr=%.4f, alpha=%.2f, beta=%.2f\n', Fr_to_pass, alpha_to_pass, beta_to_pass);

%% CALL process_eigenfunction

Xi_grid = []; B_grid = []; Psi_grid = []; 
try
    % Corrected call: Pass the full eigenvector matrix and the selected index
    [Xi_grid, B_grid, Psi_grid] = process_eigenfunction(...
        V_full_converged, eigenfunction_index_to_plot, ... % Changed here
        Kx_mesh, Kz_mesh, ...
        Fr_to_pass, alpha_to_pass, beta_to_pass, ...
        Lx_computational, Lz_computational);
    fprintf('Eigenfunction processed successfully.\n');
catch ME_process
    error('Error calling process_eigenfunction: %s\nCheck if process_eigenfunction.m and its dependencies are in ''%s'' and error-free. Also check argument order.', ...
        ME_process.message, parentDir);
end


%% PREPARE GRID FOR PLOTTING
xp_edges = linspace(0, Lx_computational, Nx + 1);
zp_edges = linspace(0, Lz_computational, Nz + 1);

%% CREATE PLOTS
fprintf('Generating plots for Xi, B, Psi...\n');
fig = figure('WindowState', 'maximized', 'Color', 'w');
cmap = slanCM('magma');

% Plot Xi (Vorticity)
subplot(3,1,1);
pcolor(xp_edges/Fr, zp_edges, [real(Xi_grid), real(Xi_grid(:,1)); real(Xi_grid(Nz,:)), real(Xi_grid(Nz,1))]); 
colormap(cmap)
shading interp;
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\xi$';
cb.Label.FontSize = plot_fontsize;
cb.FontSize = plot_fontsize - 4; 
cb.LineWidth = plot_borderwidth / 2; 
cb.TickLabelInterpreter = 'latex';
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', plot_fontsize-2);
ax1 = gca;
ax1.FontSize = plot_fontsize - 4;
ax1.LineWidth = plot_borderwidth/2; 
ax1.TickLabelInterpreter = 'latex';
daspect([1 1 1])

% Plot B (Buoyancy)
subplot(3,1,2);
pcolor(xp_edges/Fr, zp_edges, [real(B_grid), real(B_grid(:,1)); real(B_grid(Nz,:)), real(B_grid(Nz,1))]); 
colormap(cmap)
shading interp;
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$b$';
cb.Label.FontSize = plot_fontsize;
cb.FontSize = plot_fontsize - 4; 
cb.LineWidth = plot_borderwidth / 2; 
cb.TickLabelInterpreter = 'latex';
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', plot_fontsize-2);
ax1 = gca;
ax1.FontSize = plot_fontsize - 4;
ax1.LineWidth = plot_borderwidth/2; 
ax1.TickLabelInterpreter = 'latex';
daspect([1 1 1])

% Plot Psi (Streamfunction)
subplot(3,1,3);
pcolor(xp_edges/Fr, zp_edges, [real(Psi_grid), real(Psi_grid(:,1)); real(Psi_grid(Nz,:)), real(Psi_grid(Nz,1))]); 
colormap(cmap)
shading interp;
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\psi$';
cb.Label.FontSize = plot_fontsize;
cb.FontSize = plot_fontsize - 4; 
cb.LineWidth = plot_borderwidth / 2; 
cb.TickLabelInterpreter = 'latex';
xlabel('$x/\mathrm{Fr}$', 'Interpreter', 'latex', 'FontSize', plot_fontsize-2);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', plot_fontsize-2);
ax1 = gca;
ax1.FontSize = plot_fontsize - 4;
ax1.LineWidth = plot_borderwidth/2; 
ax1.TickLabelInterpreter = 'latex';
daspect([1 1 1])


%% SAVE PLOT (OPTIONAL)
outputPlotFileName = sprintf('Eigenfunction_Reb%s_Alpha%s_Index%d.pdf', ...
                             Reb_str_for_name, alpha_str_for_name, eigenfunction_index_to_plot);
outputPlotFilePath = fullfile(scriptDir, outputPlotFileName); 

try
    fprintf('Saving plot to %s...\n', outputPlotFilePath);
    if exist('exportgraphics', 'file')
        exportgraphics(fig, outputPlotFilePath, 'ContentType', 'vector', 'Resolution', 300);
    else
        print_filename = outputPlotFilePath;
        if ~endsWith(print_filename, '.pdf')
            [~, name_pf, ~] = fileparts(print_filename);
            print_filename = fullfile(scriptDir, [name_pf '.pdf']);
        end
        print(fig, print_filename, '-dpdf', '-r300', '-vector'); 
    end
    fprintf('Plot saved successfully: %s\n', outputPlotFilePath);
catch ME_save
    fprintf('Error saving plot: %s\n', ME_save.message);
end

fprintf('\nFinished script: plot_selected_eigenfunction_v2.m\n');