% SCRIPT: generate_residual_report_excel.m
%
% This script trawls solution directories, loads eigenvalue and residual data,
% and generates a single Excel (.xlsx) file. Each unique Reb value will
% have its data presented on a separate sheet within the workbook.
%
% Assumes the eigenvalues in the .mat files were generated using Fr = 0.01.
%
% Expected directory structure:
% <project_root_directory>/
%   secondary_stability_solutions_maxit10000/
%     RebX.XX_alphaY.YY/
%       spectrum_RebX.XX_alphaY.YY.mat (must contain 'eigvals' and 'residuals_norm')
%     ...
%   postprocessing/
%     generate_residual_report_excel.m (this script)
%

clear;
clc;
fprintf('Starting Excel report generation script...\n');

% --- Configuration ---
% This script assumes that the Fr used to generate the eigvals in the .mat files
% was 0.01. The table will display these eigenvalues directly.
Fr_assumed_for_loaded_eigvals = 0.01; 

fprintf('Assuming loaded eigenvalues are already scaled by Fr = %f.\n', Fr_assumed_for_loaded_eigvals);

% Determine paths relative to this script's location
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir) 
    scriptDir = pwd; 
    fprintf('Assuming script is running from: %s (postprocessing folder)\n', scriptDir);
end
solutionsParentDir = fullfile(scriptDir, '..', 'secondary_stability_solutions_maxit10000');
outputExcelDir = scriptDir; % Save .xlsx file in the current (postprocessing) folder
outputExcelFileName = 'All_Reb_Eigenvalue_Summary.xlsx';
excelFullFilePath = fullfile(outputExcelDir, outputExcelFileName);


if ~isfolder(outputExcelDir)
    mkdir(outputExcelDir);
    fprintf('Created output Excel directory: %s\n', outputExcelDir);
end

% --- Find and Parse Solution Files ---
if ~isfolder(solutionsParentDir)
    error(['Base solutions directory not found: %s\n' ...
           'Please ensure the script is in the ''postprocessing'' folder ' ...
           'and ''secondary_stability_solutions_maxit10000'' is a sibling folder.'], solutionsParentDir);
end

listing = dir(solutionsParentDir);
subDirs = {listing([listing.isdir]).name};
subDirs = subDirs(~ismember(subDirs,{'.','..'})); 

if isempty(subDirs)
    fprintf('No subdirectories found in %s.\n', solutionsParentDir);
    return;
end

fprintf('Found %d potential solution subdirectories.\n', length(subDirs));

dirPattern = '^Reb(\d{1,3}\.\d{2})_alpha(\d+\.\d{2})$';
allParsedData = struct('Reb', {}, 'RebStr', {}, 'alpha', {}, 'alphaStr', {}, ...
                       'eigvals_loaded', {}, 'residuals_norm', {});
fileCounter = 0;

for i = 1:length(subDirs)
    dirName = subDirs{i};
    tokens = regexp(dirName, dirPattern, 'tokens');
    if isempty(tokens)
        continue;
    end

    Reb_str = tokens{1}{1};
    alpha_str = tokens{1}{2};
    Reb_val = str2double(Reb_str);
    alpha_val = str2double(alpha_str);

    matFileName = sprintf('spectrum_%s.mat', dirName);
    matFilePath = fullfile(solutionsParentDir, dirName, matFileName);

    if ~isfile(matFilePath)
        continue;
    end

    try
        loadedData = load(matFilePath, 'eigvals', 'residuals_norm');
        if ~isfield(loadedData, 'eigvals') || ~isfield(loadedData, 'residuals_norm')
            fprintf('  Warning: MAT file %s missing ''eigvals'' or ''residuals_norm''. Skipping.\n', matFileName);
            continue;
        end
        if isempty(loadedData.eigvals) || isempty(loadedData.residuals_norm)
             fprintf('  Warning: MAT file %s has empty ''eigvals'' or ''residuals_norm''. Skipping.\n', matFileName);
            continue;
        end
        if length(loadedData.eigvals) ~= length(loadedData.residuals_norm)
            fprintf('  Warning: MAT file %s has mismatched lengths for ''eigvals'' and ''residuals_norm''. Skipping.\n', matFileName);
            continue;
        end
        
        fileCounter = fileCounter + 1;
        allParsedData(fileCounter).Reb = Reb_val;
        allParsedData(fileCounter).RebStr = Reb_str; 
        allParsedData(fileCounter).alpha = alpha_val;
        allParsedData(fileCounter).alphaStr = alpha_str; 
        allParsedData(fileCounter).eigvals_loaded = loadedData.eigvals(:); 
        allParsedData(fileCounter).residuals_norm = loadedData.residuals_norm(:); 
    catch ME
        fprintf('  Error loading or parsing MAT file %s: %s. Skipping.\n', matFilePath, ME.message);
    end
end

if isempty(allParsedData)
    fprintf('No valid solution files found to process.\n');
    return;
end

fprintf('\nSuccessfully parsed data from %d files.\n', fileCounter);

% Delete old Excel file to start fresh
if isfile(excelFullFilePath)
    delete(excelFullFilePath);
    fprintf('Deleted existing Excel file: %s\n', excelFullFilePath);
end

% --- Group Data by Reb and Write to Excel ---
uniqueRebs = unique([allParsedData.Reb]);
uniqueRebs = sort(uniqueRebs); 
fprintf('Found %d unique Reb values to generate Excel sheets for.\n', length(uniqueRebs));

for iReb = 1:length(uniqueRebs)
    currentReb = uniqueRebs(iReb);
    currentRebStr = sprintf('%.2f', currentReb); 
    
    % Create a valid sheet name (max 31 chars, no invalid chars)
    % Replacing '.' with 'p' for Reb in sheet name is a common practice.
    sheetName = sprintf('Reb_%s', strrep(currentRebStr, '.', 'p'));
    if length(sheetName) > 31
        sheetName = sheetName(1:31); % Truncate if too long
    end
    
    fprintf('\nGenerating sheet "%s" for Reb = %s...\n', sheetName, currentRebStr);

    dataForThisReb = allParsedData([allParsedData.Reb] == currentReb);
    [~, sortIdx] = sort([dataForThisReb.alpha]); % Sort by alpha
    dataForThisReb = dataForThisReb(sortIdx);

    % Prepare data for the table
    tableHeaders = {'alpha', 'Re_lambda_Fr', 'Im_lambda_Fr', 'Residual_Norm'};
    tableContent = []; % Initialize as an empty matrix or cell array

    for k_alpha = 1:length(dataForThisReb) 
        alphaEntry = dataForThisReb(k_alpha);
        
        % Eigenvalues are assumed to be already scaled by Fr_assumed_for_loaded_eigvals (0.01)
        eigvals_for_table = alphaEntry.eigvals_loaded;
        residuals_for_table = alphaEntry.residuals_norm;
        
        numEigsForAlpha = length(eigvals_for_table);
        
        alphaCol = repmat(alphaEntry.alpha, numEigsForAlpha, 1);
        reLambdaCol = real(eigvals_for_table);
        imLambdaCol = imag(eigvals_for_table);
        
        currentAlphaData = [alphaCol, reLambdaCol, imLambdaCol, residuals_for_table];
        if isempty(tableContent)
            tableContent = currentAlphaData;
        else
            tableContent = [tableContent; currentAlphaData]; %#ok<AGROW>
        end
    end

    if isempty(tableContent)
        fprintf('  No data to write to Excel for Reb = %s (Sheet "%s"). Skipping this sheet.\n', currentRebStr, sheetName);
        continue;
    end

    % Convert to MATLAB table
    T = array2table(tableContent, 'VariableNames', tableHeaders);

    try
        writetable(T, excelFullFilePath, 'Sheet', sheetName);
        fprintf('  Successfully written data for Reb %s to sheet "%s" in %s\n', ...
                currentRebStr, sheetName, outputExcelFileName);
    catch ME_excel
        fprintf('  Error writing to Excel sheet "%s" for Reb = %s: %s\n', ...
                sheetName, currentRebStr, ME_excel.message);
        fprintf('  Make sure the Excel file is not open or locked by another process.\n');
    end
end

fprintf('\n--- Excel Report Generation Complete ---\n');
if exist(excelFullFilePath, 'file')
    fprintf('Excel file created: %s\n', excelFullFilePath);
else
    fprintf('Excel file was not created (perhaps no data was processed).\n');
end