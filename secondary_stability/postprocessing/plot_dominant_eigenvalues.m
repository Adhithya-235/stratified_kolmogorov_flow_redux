% SCRIPT: plot_dominant_eigenvalues_by_Reb.m
%
% This script trawls solution folders, loads converged eigenpairs,
% filters for Re(lambda) >= 0, identifies the dominant eigenvalue (largest 
% real part) for each case (Reb, alpha pair). It then generates a scatter 
% plot of these dominant eigenvalues in the complex plane. 
% - For each alpha, lines connect dominant eigenvalues from different Rebs.
% - Markers are colored by their Reb value.
% - Complex conjugates of dominant eigenvalues (if Im~=0) are also plotted.
%
% Assumes eigenvalues in loaded .mat files are scaled by Fr = 0.01.
%

clear;
clc;
close all; 

%% SCRIPT INITIALIZATION & CONFIGURATION

addpath('../../utility_belt');
fprintf('Starting script: plot_dominant_eigenvalues_by_Reb.m\n');

% Configuration for plotting
Fr_assumed_for_eigenvalues = 0.01; 
plot_fontsize = 20;
plot_borderwidth = 3;
marker_size = 80; 
colormap_to_use = slanCM('deep'); 
output_plot_filename = 'Dominant_Eigenvalues.pdf'; 

% Define fixed axis limits (adjust as needed)
x_plot_limits = [0, 0.15]; % Example: Focus near origin and positive real
y_plot_limits = [-1.5, 1.5]; % Example

%% PATH AND DIRECTORY SETUP
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir) 
    scriptDir = pwd; 
    fprintf('Assuming script is running from: %s (postprocessing folder)\n', scriptDir);
end
solutionsParentDir = fullfile(scriptDir, '..', 'secondary_stability_solutions_maxit10000');
outputPlotDir = scriptDir; 

if ~isfolder(outputPlotDir)
    mkdir(outputPlotDir);
    fprintf('Created output plot directory: %s\n', outputPlotDir);
end

if ~isfolder(solutionsParentDir)
    error(['Base solutions directory not found: %s\n' ...
           'Ensure script is in ''postprocessing'' and solutions folder is sibling.'], solutionsParentDir);
end

%% LOCATE AND PARSE CONVERGED DATA FILES
listing = dir(solutionsParentDir);
subDirs = {listing([listing.isdir]).name};
subDirs = subDirs(~ismember(subDirs,{'.','..'})); 

if isempty(subDirs)
    fprintf('No subdirectories found in %s.\n', solutionsParentDir);
    return;
end
fprintf('Found %d potential solution subdirectories to check.\n', length(subDirs));

dirPattern = '^Reb(\d{1,3}\.\d{2})_alpha(\d+\.\d{2})$'; 

% Store all relevant data points: Reb, alpha, dominant_eig_with_Re_ge_0
all_points_data = struct('Reb', [], 'alpha', [], 'dominant_eig', []);
data_idx = 0;

for i = 1:length(subDirs)
    dirName = subDirs{i}; 
    tokens = regexp(dirName, dirPattern, 'tokens');
    if isempty(tokens)
        continue; 
    end

    Reb_str = tokens{1}{1};      
    alpha_str = tokens{1}{2};    
    Reb_val_current = str2double(Reb_str); 
    alpha_val_current = str2double(alpha_str);

    Reb_str_for_name = strrep(Reb_str, '.', 'p');     
    alpha_str_for_name = strrep(alpha_str, '.', 'p'); 
    
    convergedMatFilename = sprintf('converged_Reb%s_alpha%s.mat', Reb_str_for_name, alpha_str_for_name);
    convergedMatFilePath = fullfile(solutionsParentDir, dirName, convergedMatFilename);
    dynamicVarName = sprintf('converged_Reb%s_alpha%s', Reb_str_for_name, alpha_str_for_name);

    if isfile(convergedMatFilePath)
        fprintf('  Processing file: %s\n', convergedMatFilePath);
        try
            loadedFile = load(convergedMatFilePath);
            if isfield(loadedFile, dynamicVarName)
                dataStruct = loadedFile.(dynamicVarName);
                
                if isfield(dataStruct, 'converged_eigenvalues') && ~isempty(dataStruct.converged_eigenvalues)
                    % Filter for eigenvalues with non-negative real part
                    eigs_to_consider = dataStruct.converged_eigenvalues(real(dataStruct.converged_eigenvalues) >= -1e-12); % allow very small negative due to precision

                    if ~isempty(eigs_to_consider)
                        % Find the dominant eigenvalue (largest real part among those with Re >= 0)
                        [~, idx_dominant] = max(real(eigs_to_consider));
                        dominant_eig = eigs_to_consider(idx_dominant);
                        
                        data_idx = data_idx + 1;
                        all_points_data(data_idx).Reb = dataStruct.Reb;
                        all_points_data(data_idx).alpha = dataStruct.alpha;
                        all_points_data(data_idx).dominant_eig = dominant_eig;
                    else
                        fprintf('    No converged eigenvalues with Re >= 0 found for Reb %s, alpha %s.\n', Reb_str, alpha_str);
                    end
                else
                    fprintf('    ''converged_eigenvalues'' field empty or not found in struct ''%s''.\n', dynamicVarName);
                end
            else
                fprintf('    Dynamically named struct ''%s'' not found in MAT file. Skipping.\n', dynamicVarName);
            end
        catch ME
            fprintf('    Error loading or processing MAT file %s: %s. Skipping.\n', convergedMatFilePath, ME.message);
        end
    end
end

if data_idx == 0
    fprintf('\nNo dominant eigenvalues (with Re>=0) found from any processed files. Cannot generate plot.\n');
    return;
end

fprintf('\nExtracted %d dominant eigenvalues (with Re>=0) for plotting.\n', data_idx);

%% PREPARE DATA FOR SCATTER PLOT (MARKERS AND CONJUGATES)
scatter_eigs_real = [];
scatter_eigs_imag = [];
scatter_rebs_color = [];
imag_tolerance = 1e-9; 

for k_point = 1:length(all_points_data)
    eig_val = all_points_data(k_point).dominant_eig;
    reb_val = all_points_data(k_point).Reb;
    
    scatter_eigs_real = [scatter_eigs_real; real(eig_val)]; %#ok<AGROW>
    scatter_eigs_imag = [scatter_eigs_imag; imag(eig_val)]; %#ok<AGROW>
    scatter_rebs_color = [scatter_rebs_color; reb_val]; %#ok<AGROW>
    
    if abs(imag(eig_val)) > imag_tolerance
        scatter_eigs_real = [scatter_eigs_real; real(eig_val)]; %#ok<AGROW>
        scatter_eigs_imag = [scatter_eigs_imag; -imag(eig_val)];%#ok<AGROW>
        scatter_rebs_color = [scatter_rebs_color; reb_val]; %#ok<AGROW>
    end
end

%% CREATE PLOT
fprintf('Generating plot...\n');
fig = figure('WindowState', 'maximized', 'Color', 'w');
ax = gca; % Get current axes
hold(ax, 'on');

%% SCATTER PLOT FOR DOMINANT EIGENVALUES AND THEIR CONJUGATES (MARKERS)
scatter(ax, scatter_eigs_real, scatter_eigs_imag, marker_size, scatter_rebs_color, ...
    'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.75);

%% COLORMAP AND COLORBAR (FOR REB)
colormap(ax, colormap_to_use); 
cb = colorbar(ax);
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\mathrm{Re_b}$';
cb.Label.FontSize = plot_fontsize;
cb.FontSize = plot_fontsize - 4; 
cb.LineWidth = plot_borderwidth / 2; 
cb.TickLabelInterpreter = 'latex';

%% AXES LABELS AND TITLE
xlabel_str = '$\mathrm{Fr}\,\mathcal{R}\left[\sigma\right]$';
ylabel_str = '$\mathrm{Fr}\,\mathcal{I}\left[\sigma\right]$';
xlabel(xlabel_str, 'Interpreter', 'latex', 'FontSize', plot_fontsize);
ylabel(ylabel_str, 'Interpreter', 'latex', 'FontSize', plot_fontsize);

%% AXES PROPERTIES
ax.FontSize = plot_fontsize - 2; 
ax.LineWidth = plot_borderwidth;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = ':';
ax.TickLabelInterpreter = 'latex'; % Already set for colorbar, usually applies to axes too

axis(ax, 'square'); % Makes the plot box square

% Apply fixed limits
xlim(ax, x_plot_limits);
ylim(ax, y_plot_limits);
hold(ax, 'off');

%% SAVE PLOT
outputPlotFilePath = fullfile(outputPlotDir, output_plot_filename);
try
    fprintf('Saving plot to %s...\n', outputPlotFilePath);
    if exist('exportgraphics', 'file')
        exportgraphics(fig, outputPlotFilePath, 'ContentType', 'vector', 'Resolution', 300);
    else
        print_filename = outputPlotFilePath;
        if ~endsWith(print_filename, '.pdf') % print needs .pdf for -dpdf
            [~, name, ~] = fileparts(print_filename);
            print_filename = fullfile(outputPlotDir, [name '.pdf']);
            fprintf('  Saving as PDF for print command: %s\n', print_filename);
        end
        print(fig, print_filename, '-dpdf', '-r300', '-vector'); 
    end
    fprintf('Plot saved successfully.\n');
catch ME_save
    fprintf('Error saving plot: %s\n', ME_save.message);
end

%% FINAL SUMMARY
fprintf('Finished script: plot_dominant_eigenvalues_by_Reb.m\n');