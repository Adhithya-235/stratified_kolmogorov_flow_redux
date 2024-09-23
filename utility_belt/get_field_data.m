function [t, u, w, b, p, vort, nf] = get_field_data(folder_name, data_folder, file_name, stride, svec, wrap)

% This function reads HDF5 data produced by dedalus and extracts temporal
% grid as well as 3-dimensional (x, z, t) primitive variable fields. 
% Specify date-based folder name, data file name as well as the desired 
% data series numbers to be plotted. Also, wrap=1 means the last gridpoint 
% in periodic grids is kept. Useful for integration purposes. Stride spec-
% -ifies the index interval in the final dimension of the hdf5 data.


%% FILENAME

fname = string.empty;
maxs = length(svec);
for s = 1:maxs
    fname(s) = sprintf('../%s/%s/%s/%s_s%d.h5', folder_name, data_folder, file_name, file_name, svec(s));
end

%% GET DATA FROM FILE

u    = [];
w    = [];
b    = [];
p    = [];
vort = [];
t    = [];

for s = 1:maxs
   fprintf('Reading file %d of %d.\n',svec(s),svec(maxs)) 
   t    = [t; h5read(fname(s),'/scales/sim_time',1,Inf,stride)];
   u    = cat(3, u, h5read(fname(s), '/tasks/u',[1,1,1],[Inf,Inf,Inf],[1,1,stride]));
   w    = cat(3, w, h5read(fname(s), '/tasks/w',[1,1,1],[Inf,Inf,Inf],[1,1,stride]));
   b    = cat(3, b, h5read(fname(s), '/tasks/b',[1,1,1],[Inf,Inf,Inf],[1,1,stride]));
   p    = cat(3, p, h5read(fname(s), '/tasks/p',[1,1,1],[Inf,Inf,Inf],[1,1,stride]));
   vort = cat(3, vort, h5read(fname(s), '/tasks/o',[1,1,1],[Inf,Inf,Inf],[1,1,stride]));
 end


%% WRAP

if wrap == 1
    u    = cat(2, u, u(:, 1, :));
    w    = cat(2, w, w(:, 1, :));
    b    = cat(2, b, b(:, 1, :));
    p    = cat(2, p, p(:, 1, :));
    vort = cat(2, vort, vort(:, 1, :));
end

%% DETERMINE TIMESERIES LENGTH

nf = length(t);

end