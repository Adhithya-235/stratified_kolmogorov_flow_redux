function [x, z, X, Z] = get_space_data(folder_name, data_folder, file_name, wrap)

% This function reads HDF5 data produced by dedalus and extracts the spatial grid. 
% Specify date-based folder name and data file name. Also, wrap=1 means the last 
% gridpoint in periodic grids is kept. Useful for integration purposes.

%% FILENAME 

fname = sprintf('../%s/%s/%s/%s_s1.h5', folder_name, data_folder, file_name, file_name);

%% READ x, y AND z DATA

x = h5read(fname,'/scales/x/1.0');
z = h5read(fname,'/scales/z/1.0');

%% WRAP 

if wrap == 1
    x = cat(1, x, 2*x(end)-x(end-1));
    z = cat(1, z, 2*z(end)-z(end-1));
end

%% MESHGRID

[X, Z] = meshgrid(x, z);

end