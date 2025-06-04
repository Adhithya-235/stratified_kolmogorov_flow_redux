% SCRIPT: isolate_converged_eigenpairs.m
%
% This script trawls solution directories, isolates eigenpairs with
% residual norm below a specified threshold, sorts them by the descending
% magnitude of the real part of their eigenvalues, and stores them along
% with relevant metadata in a new, dynamically named .mat file within
% each solution folder. The data within this .mat file is stored in a 
% dynamically named structure variable (e.g., converged_Reb1p23_alpha0p45).
%
% Assumes the eigenvalues in the source .mat files were generated 
% using Fr = 0.01.
%
% Expected directory structure:
% <project_root_directory>/
%   secondary_stability_solutions_maxit10000/
%     RebX.XX_alphaY.YY/
%       spectrum_RebX.XX_alphaY.YY.mat (source file)
%       converged_RebXpXX_alphaYpYY.mat (will be created by this script,
%                                       containing one dynamically named 
%                                       structure variable)
%     ...
%   postprocessing/
%     isolate_converged_eigenpairs.m (this script)
%

%% SCRIPT INITIALIZATION
clear;
clc;
fprintf('Starting script: isolate_converged_eigenpairs.m (v2 - dynamic filename)...\n');

%% USER CONFIGURATION
% This Fr value is assumed for the scaling of 'eigvals' loaded from source files.
% It will also be recorded in the output metadata.
Fr_scaling_of_eigenvalues = 0.01; 
residual_threshold = 5e-8; % Residual norm threshold

fprintf('Using Fr scaling for eigenvalues: %f\n', Fr_scaling_of_eigenvalues);
fprintf('Residual norm threshold for convergence: %e\n', residual_threshold);

%% PATH AND DIRECTORY SETUP
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir) 
    scriptDir = pwd; 
    fprintf('Assuming script is running from: %s (postprocessing folder)\n', scriptDir);
end
solutionsParentDir = fullfile(scriptDir, '..', 'secondary_stability_solutions_maxit10000');

if ~isfolder(solutionsParentDir)
    error(['Base solutions directory not found: %s\n' ...
           'Please ensure the script is in the ''postprocessing'' folder ' ...
           'and ''secondary_stability_solutions_maxit10000'' is a sibling folder.'], solutionsParentDir);
end

%% LOCATE SOLUTION SUBDIRECTORIES
listing = dir(solutionsParentDir);
subDirs = {listing([listing.isdir]).name};
subDirs = subDirs(~ismember(subDirs,{'.','..'})); 

if isempty(subDirs)
    fprintf('No subdirectories found in %s.\n', solutionsParentDir);
    return;
end
fprintf('Found %d potential solution subdirectories to process.\n', length(subDirs));

%% INITIALIZE COUNTERS AND PATTERNS
dirPattern = '^Reb(\d{1,3}\.\d{2})_alpha(\d+\.\d{2})$'; 
processedCount = 0;
skippedCount = 0;
noConvergedFoundCount = 0;

%% MAIN PROCESSING LOOP (PER SOLUTION FOLDER)
for i = 1:length(subDirs)
    dirName = subDirs{i}; 
    currentSolutionFolderPath = fullfile(solutionsParentDir, dirName);
    
    fprintf('\n%% PROCESSING DIRECTORY: %s\n', dirName);

    tokens = regexp(dirName, dirPattern, 'tokens');
    if isempty(tokens)
        fprintf('  Skipping directory (pattern mismatch): %s\n', dirName);
        skippedCount = skippedCount + 1;
        continue;
    end

    Reb_str = tokens{1}{1};      
    alpha_str = tokens{1}{2};    
    Reb_val = str2double(Reb_str);
    alpha_val = str2double(alpha_str);

    Reb_str_for_name = strrep(Reb_str, '.', 'p');     
    alpha_str_for_name = strrep(alpha_str, '.', 'p'); 

    sourceMatFileName = sprintf('spectrum_%s.mat', dirName);
    sourceMatFilePath = fullfile(currentSolutionFolderPath, sourceMatFileName);

    if ~isfile(sourceMatFilePath)
        fprintf('  Skipping (Source MAT file not found): %s\n', sourceMatFilePath);
        skippedCount = skippedCount + 1;
        continue;
    end

    %% LOAD SOURCE DATA
    fprintf('  Loading source MAT file: %s\n', sourceMatFileName);
    loadedData = []; 
    try
        loadedData = load(sourceMatFilePath, 'eigvals', 'eigvecs', 'residuals_norm', 'meta');
    catch ME
        fprintf('  Error loading source MAT file %s: %s. Skipping.\n', sourceMatFileName, ME.message);
        skippedCount = skippedCount + 1;
        continue;
    end

    %% VALIDATE LOADED DATA
    requiredFields = {'eigvals', 'eigvecs', 'residuals_norm', 'meta'};
    missingFields = setdiff(requiredFields, fieldnames(loadedData));
    if ~isempty(missingFields)
        fprintf('  Skipping MAT file (missing top-level fields: %s).\n', strjoin(missingFields, ', '));
        skippedCount = skippedCount + 1;
        continue;
    end
    if ~isfield(loadedData.meta, 'resolution') || ~isfield(loadedData.meta, 'domainsize')
        fprintf('  Skipping MAT file (meta struct missing ''resolution'' or ''domainsize'').\n');
        skippedCount = skippedCount + 1;
        continue;
    end
    if isempty(loadedData.eigvals) || isempty(loadedData.eigvecs) || isempty(loadedData.residuals_norm)
        fprintf('  Skipping MAT file (''eigvals'', ''eigvecs'', or ''residuals_norm'' is empty).\n');
        skippedCount = skippedCount + 1;
        continue;
    end
    if length(loadedData.eigvals) ~= length(loadedData.residuals_norm) || ...
       size(loadedData.eigvecs, 2) ~= length(loadedData.eigvals)
        fprintf('  Skipping MAT file (data dimensions mismatch).\n');
        skippedCount = skippedCount + 1;
        continue;
    end

    %% FILTER CONVERGED EIGENPAIRS
    converged_indices = loadedData.residuals_norm < residual_threshold;

    if ~any(converged_indices)
        fprintf('  No eigenvalues found below residual threshold %e for %s.\n', residual_threshold, dirName);
        noConvergedFoundCount = noConvergedFoundCount + 1;
        continue; 
    end

    eigvals_converged_raw = loadedData.eigvals(converged_indices);
    eigvecs_converged_raw = loadedData.eigvecs(:, converged_indices);
    resnorms_converged_raw = loadedData.residuals_norm(converged_indices);
    
    fprintf('  Found %d eigenpairs below residual threshold.\n', length(eigvals_converged_raw));

    %% SORT CONVERGED EIGENPAIRS
    [~, sort_order] = sort(real(eigvals_converged_raw), 'descend');
    
    final_eigvals = eigvals_converged_raw(sort_order);
    final_eigvecs = eigvecs_converged_raw(:, sort_order);
    final_resnorms = resnorms_converged_raw(sort_order);

    %% PREPARE AND SAVE CONVERGED DATA
    structureToSave = struct(...
        'Reb', Reb_val, ...
        'alpha', alpha_val, ...
        'original_resolution', loadedData.meta.resolution, ...
        'original_domainsize', loadedData.meta.domainsize, ...
        'Fr_scaling_of_eigenvalues', Fr_scaling_of_eigenvalues, ...
        'residual_norm_threshold_used', residual_threshold, ...
        'converged_eigenvalues', final_eigvals, ...
        'converged_eigenvectors', final_eigvecs, ...
        'converged_residual_norms', final_resnorms ...
    );

    dynamicVarName = sprintf('converged_Reb%s_alpha%s', Reb_str_for_name, alpha_str_for_name);
    eval([dynamicVarName ' = structureToSave;']);

    % Define dynamic output MAT file name
    outputMatFileName = sprintf('%s.mat', dynamicVarName); 
    outputMatFilePath = fullfile(currentSolutionFolderPath, outputMatFileName);
    
    try
        save(outputMatFilePath, dynamicVarName); % Save only the dynamically named variable
        fprintf('  Saved %d converged eigenpairs to %s (variable: %s)\n', ...
                length(final_eigvals), outputMatFileName, dynamicVarName);
        processedCount = processedCount + 1;
    catch ME_save
        fprintf('  Error saving converged data to %s: %s\n', outputMatFileName, ME_save.message);
        skippedCount = skippedCount + 1;
    end
    
    clear(dynamicVarName); % Clean up from workspace

end

%% FINAL SUMMARY
fprintf('\n%% SCRIPT EXECUTION SUMMARY\n');
fprintf('Successfully processed %d solution folders and saved converged data.\n', processedCount);
fprintf('Skipped %d folders due to missing data, errors, or unmet criteria.\n', skippedCount);
fprintf('%d folders had no eigenpairs meeting the convergence threshold.\n', noConvergedFoundCount);
fprintf('Finished script: isolate_converged_eigenpairs.m\n');