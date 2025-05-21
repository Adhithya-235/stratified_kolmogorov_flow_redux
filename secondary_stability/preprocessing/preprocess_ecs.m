clear 
close all
clc

%% Directory containing .mat files

data_dir = '../import_test';
files = dir(fullfile(data_dir, '*.mat'));

for file = files'
    
    filepath = fullfile(data_dir, file.name);
    data = load(filepath);
    modified = false;

    %% --- Step 0: Check for Re_b ---
    
    if ~isfield(data, 'Re_b')
        if isfield(data, 'Buoyancy') && isfield(data, 'z')
            
            %% Make sure z is column vector for proper broadcasting
            
            z_vec = data.z(:);
            
            %% Subtract z from each column of Buoyancy
            
            data.Buoyancy = data.Buoyancy - z_vec;
            modified = true;
        else
            
            warning('Missing Buoyancy or z in file %s. Skipping z subtraction.', file.name);
        
        end
    end

    %% --- Step 1 & 2: Check and interpolate arrays ---
    
    fields_2D = {'Buoyancy', 'omega', 'psi'};
    fields_1D = {'x', 'z'};
    
    for f = fields_2D
        fname = f{1};
        if isfield(data, fname)
            A = data.(fname);
            if isequal(size(A), [256, 256])
                data.(fname) = interpft(interpft(A, 128, 1), 128, 2);
                modified = true;
            end
        end
    end

    for f = fields_1D
        fname = f{1};
        if isfield(data, fname)
            A = data.(fname);
            if isvector(A) && length(A) == 256
                data.(fname) = interpft(A, 128);
                modified = true;
            end
        end
    end

    %% --- Save back if anything changed ---
    
    if modified
        save(filepath, '-struct', 'data');
        fprintf('Modified and saved: %s\n', file.name);
    else
        fprintf('No changes needed: %s\n', file.name);
    end
end
