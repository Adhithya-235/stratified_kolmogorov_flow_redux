% This script trawls solution folders, loads the MAT-files containing
% converged eigenpairs (created by isolate_converged_eigenpairs.m). [cite: 1]
% For each alpha, it plots the real parts of the first five converged 
% eigenvalues against Reb. [cite: 2]
% Each of the five eigenvalues forms a separate line. [cite: 3]
%
% Assumes eigenvalues in loaded .mat files are scaled by Fr = 0.01.
% This version creates a separate PDF file for each unique alpha value.

clear;
clc;
close all; 

%% SCRIPT INITIALIZATION & CONFIGURATION

addpath('../../utility_belt'); % User's added path
fprintf('Starting script: plot_growth_rates_per_alpha.m\n'); % Updated script name in log

% Configuration
Fr_assumed_for_eigenvalues = 0.01; 
num_eigenvalues_to_plot = 3; % Plot up to the first 5 [cite: 5]
plot_fontsize = 20;
plot_borderwidth = 3;
plot_linewidth = 3;
marker_symbol = 'o';
% output_plot_filename = 'Growth_Rates.pdf'; % Old single filename[cite: 6], now dynamic per alpha

% Define fixed Y-axis limits for Re(lambda) (adjust as needed for your data)
y_plot_limits_re_lambda = [-0.15, 0.15]; % Example limits [cite: 7]

%% PATH AND DIRECTORY SETUP
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir) 
    scriptDir = pwd; 
    fprintf('Assuming script is running from: %s (postprocessing folder)\n', scriptDir); % [cite: 8]
end
solutionsParentDir = fullfile(scriptDir, '..', 'secondary_stability_solutions_maxit10000');
outputPlotDir = scriptDir; % Save PDFs in the current (postprocessing) folder
if ~isfolder(outputPlotDir) % [cite: 9]
    mkdir(outputPlotDir);
    fprintf('Created output plot directory: %s\n', outputPlotDir); % [cite: 9]
end

if ~isfolder(solutionsParentDir) % [cite: 10]
    error(['Base solutions directory not found: %s\n' ...
           'Ensure script is in ''postprocessing'' and solutions folder is sibling.'], solutionsParentDir); % [cite: 10]
end

%% LOCATE AND PARSE CONVERGED DATA FILES
listing = dir(solutionsParentDir);
subDirs = {listing([listing.isdir]).name};
subDirs = subDirs(~ismember(subDirs,{'.','..'})); 
if isempty(subDirs) % [cite: 12]
    fprintf('No subdirectories found in %s.\n', solutionsParentDir); % [cite: 12]
    return; % [cite: 12]
end
fprintf('Found %d potential solution subdirectories to check.\n', length(subDirs));
dirPattern = '^Reb(\d{1,3}\.\d{2})_alpha(\d+\.\d{2})$';  % [cite: 13]

% Store all relevant data points: Reb, alpha, list of top N eigenvalues
collected_data = struct('Reb', [], 'alpha', [], 'top_N_eigenvalues', {}); 
data_idx = 0; % [cite: 14]

for i = 1:length(subDirs)
    dirName = subDirs{i}; 
    tokens = regexp(dirName, dirPattern, 'tokens'); 
    if isempty(tokens) % [cite: 15]
        continue;  % [cite: 15]
    end

    Reb_str = tokens{1}{1};      
    alpha_str = tokens{1}{2};    
    Reb_val_current = str2double(Reb_str);  % [cite: 16]
    alpha_val_current = str2double(alpha_str); % [cite: 16]

    Reb_str_for_name = strrep(Reb_str, '.', 'p');     
    alpha_str_for_name = strrep(alpha_str, '.', 'p'); 
    
    convergedMatFilename = sprintf('converged_Reb%s_alpha%s.mat', Reb_str_for_name, alpha_str_for_name); % [cite: 17]
    convergedMatFilePath = fullfile(solutionsParentDir, dirName, convergedMatFilename);
    dynamicVarName = sprintf('converged_Reb%s_alpha%s', Reb_str_for_name, alpha_str_for_name);
    
    if isfile(convergedMatFilePath) % [cite: 18]
        fprintf('  Processing file: %s\n', convergedMatFilePath); % [cite: 18]
        try % [cite: 19]
            loadedFile = load(convergedMatFilePath); % [cite: 19]
            if isfield(loadedFile, dynamicVarName) % [cite: 20]
                dataStruct = loadedFile.(dynamicVarName); % [cite: 20]
                if isfield(dataStruct, 'converged_eigenvalues') && ~isempty(dataStruct.converged_eigenvalues) % [cite: 21]
                    num_available_eigs = length(dataStruct.converged_eigenvalues); % [cite: 21]
                    num_to_take = min(num_available_eigs, num_eigenvalues_to_plot); % [cite: 22]
                    
                    first_n_eigs = dataStruct.converged_eigenvalues(1:num_to_take); % [cite: 22]
                    
                    data_idx = data_idx + 1;
                    collected_data(data_idx).Reb = dataStruct.Reb;
                    collected_data(data_idx).alpha = dataStruct.alpha;
                    collected_data(data_idx).top_N_eigenvalues = first_n_eigs; % [cite: 23]
                else
                    fprintf('    ''converged_eigenvalues'' field empty or not found in struct ''%s''.\n', dynamicVarName); % [cite: 24]
                end
            else
                fprintf('    Dynamically named struct ''%s'' not found in MAT file. Skipping.\n', dynamicVarName); % [cite: 25]
            end
        catch ME
            fprintf('    Error loading or processing MAT file %s: %s. Skipping.\n', convergedMatFilePath, ME.message); % [cite: 26]
        end
    end
end

if data_idx == 0
    fprintf('\nNo converged eigenvalues found from any processed files. Cannot generate plot.\n'); % [cite: 27]
    return; % [cite: 27]
end

fprintf('\nExtracted data for %d (Reb, alpha) pairs.\n', data_idx);

%% PREPARE AND CREATE PLOTS (ONE SEPARATE PDF PER ALPHA)
unique_alphas = sort(unique([collected_data.alpha]));
fprintf('Found %d unique alpha values to generate plots for.\n', length(unique_alphas)); % [cite: 28]

% outputPlotFilePath = fullfile(outputPlotDir, output_plot_filename); % Old single file path
% if isfile(outputPlotFilePath) % [cite: 29]
% delete(outputPlotFilePath); % Delete old PDF to start fresh [cite: 29]
% fprintf('Deleted existing PDF: %s\n', outputPlotFilePath); % [cite: 29]
% end
% firstPagePrinted = false; % No longer needed for separate files

mode_colors = lines(num_eigenvalues_to_plot); % Get distinct colors for the N modes [cite: 30]

for i_alpha = 1:length(unique_alphas)
    current_alpha = unique_alphas(i_alpha);
    fprintf('\nGenerating plot for alpha = %.2f...\n', current_alpha); % [cite: 31]
    
    % Filter data for current alpha
    alpha_specific_data_indices = find([collected_data.alpha] == current_alpha);
    if isempty(alpha_specific_data_indices) % [cite: 32]
        fprintf('  No data for alpha = %.2f. Skipping page.\n', current_alpha); % [cite: 32]
        continue; % [cite: 33]
    end
    current_alpha_data = collected_data(alpha_specific_data_indices);
    
    % Sort by Reb for proper line plotting
    [sorted_Rebs, sort_idx] = sort([current_alpha_data.Reb]);
    current_alpha_data_sorted = current_alpha_data(sort_idx); % [cite: 34]
    
    if isempty(sorted_Rebs)
        fprintf('  No Reb values for alpha = %.2f after sorting. Skipping page.\n', current_alpha); % [cite: 35]
        continue; % [cite: 35]
    end

    fig = figure('WindowState', 'maximized', 'Color', 'w');
    ax = gca;
    hold(ax, 'on');
    
    legend_handles = [];
    legend_labels = {}; % [cite: 36] % Though not used if DisplayName is used in plot

    for mode_num = 1:num_eigenvalues_to_plot
        real_parts_for_mode = NaN(size(sorted_Rebs)); % Initialize with NaNs [cite: 37]
        
        for k_reb = 1:length(sorted_Rebs)
            if length(current_alpha_data_sorted(k_reb).top_N_eigenvalues) >= mode_num
                real_parts_for_mode(k_reb) = real(current_alpha_data_sorted(k_reb).top_N_eigenvalues(mode_num)); % [cite: 38]
            end
        end
        
        if any(~isnan(real_parts_for_mode))
            h_plot = plot(ax, sorted_Rebs, real_parts_for_mode, ...
                ['-' marker_symbol], ...
                'Color', mode_colors(mode_num, :), ...
                'LineWidth', plot_linewidth, ... % [cite: 39]
                'MarkerSize', 8, ...
                'MarkerFaceColor', mode_colors(mode_num, :), ...
                'MarkerEdgeColor', 'k', ...
                'DisplayName', sprintf('Mode %d', mode_num)); % [cite: 39]
            legend_handles = [legend_handles, h_plot]; %#ok<AGROW> % [cite: 40]
        end
    end
    
    hold(ax, 'off');
    
    if isempty(legend_handles) % [cite: 41]
        fprintf('  No data plotted for alpha = %.2f. Skipping PDF file.\n', current_alpha); % [cite: 41]
        if ishandle(fig); close(fig); end % [cite: 42]
        continue; % [cite: 42]
    end
    
    % Axes Labels and Title
    xlabel_str = '$\mathrm{Re_b}$'; % [cite: 43]
    ylabel_str = '$\mathrm{Fr}\,\mathcal{R}\left[\sigma\right]$'; % [cite: 43]
    % Define title_str for each plot
    title_str = ['Growth Rates ($\alpha = ', num2str(current_alpha), '$)'];
    
    xlabel(ax, xlabel_str, 'Interpreter', 'latex', 'FontSize', plot_fontsize);
    ylabel(ax, ylabel_str, 'Interpreter', 'latex', 'FontSize', plot_fontsize);
    title(ax, title_str, 'Interpreter', 'latex', 'FontSize', plot_fontsize + 2); % [cite: 44] using new title_str

    % Legend
    legend(ax, legend_handles, 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', plot_fontsize-4);
    
    % Axes Properties [cite: 45]
    ax.FontSize = plot_fontsize - 2; 
    ax.LineWidth = plot_borderwidth;
    ax.Box = 'on';
    ax.XGrid = 'on'; 
    ax.YGrid = 'on'; % [cite: 46]
    ax.GridLineStyle = ':';
    ax.TickLabelInterpreter = 'latex'; 

    axis(ax, 'square'); 
    ylim(ax, y_plot_limits_re_lambda); % Apply fixed Y limits [cite: 47]
    if ~isempty(sorted_Rebs)
        % Dynamic X limits based on Reb range for this alpha [cite: 48]
        x_min_limit = min(sorted_Rebs);
        x_max_limit = max(sorted_Rebs);
        if x_min_limit == x_max_limit % Single Reb point
            x_range_padding = max(1, 0.1 * abs(x_min_limit)); % Ensure some padding
            xlim(ax, [x_min_limit - x_range_padding, x_max_limit + x_range_padding]);
        else
            padding = 0.05 * (x_max_limit - x_min_limit); % 5% padding
            if padding == 0; padding = 1; end % Ensure some padding if range is zero
            xlim(ax, [x_min_limit - padding, x_max_limit + padding]);
        end
    end

    % Define dynamic PDF filename for this alpha
    current_alpha_str_for_file = strrep(sprintf('%.2f', current_alpha), '.', 'p');
    outputPlotFileNameThisAlpha = sprintf('Growth_Rates_alpha_%s.pdf', current_alpha_str_for_file);
    outputPlotFilePathThisAlpha = fullfile(outputPlotDir, outputPlotFileNameThisAlpha);

    % Save plot to a separate PDF for this alpha
    try
        fprintf('  Creating PDF: %s\n', outputPlotFileNameThisAlpha);
        if exist('exportgraphics', 'file')
            exportgraphics(fig, outputPlotFilePathThisAlpha, 'ContentType', 'vector', 'Resolution', 300); % No 'Append'
        else
            print_filename_this_alpha = outputPlotFilePathThisAlpha;
            if ~endsWith(print_filename_this_alpha, '.pdf') % [cite: 55]
                [~, name_pf, ~] = fileparts(print_filename_this_alpha); % [cite: 55]
                print_filename_this_alpha = fullfile(outputPlotDir, [name_pf '.pdf']); % [cite: 56]
             end
            print(fig, print_filename_this_alpha, '-dpdf', '-r300', '-painters'); % No 'Append'
        end
        fprintf('    PDF for alpha %.2f saved: %s\n', current_alpha, outputPlotFileNameThisAlpha);
    catch ME_save % Changed ME_export to ME_save for clarity
        fprintf('    Error saving PDF for alpha = %.2f: %s\n', current_alpha, ME_save.message); % [cite: 53] using ME_save
        % Fallback logic from user's script was complex and part of append logic, simplified here
    end
    
    if ishandle(fig)
        close(fig); % [cite: 61]
    end
end

%% FINAL SUMMARY
fprintf('\n--- Plot Generation Complete ---');
fprintf('\nIndividual PDF files for each alpha should be available in: %s\n', outputPlotDir); % [compare with cite: 62, 63]
fprintf('Finished script: plot_growth_rates_per_alpha.m\n'); % Updated script name in log