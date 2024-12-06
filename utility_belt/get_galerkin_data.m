function [x, z, t, u, w, b, vort, nf] = get_galerkin_data(folder_name, data_folder, file_name, svec, nx, nz, modes, wrap)

% This function reads HDF5 data produced by dedalus and extracts temporal
% grid as well as 2-dimensional (z, t) galerkin-projected primitive variable fields. 
% Specify date-based folder name, data file name as well as the desired 
% data series numbers to be plotted. Also, wrap=1 means the last gridpoint 
% in periodic grids is kept. Useful for integration purposes. Stride spec-
% -ifies the index interval in the final dimension of the hdf5 data. The variable
% modes is a vector containing the mode indices used in the truncation. There are 4 variables
% produced by the simulation for each mode -- u, w, b, vorticity. Usually arranged in that order.     

%% PREPROCESSING

Froude  = 0.02;
Lx      = Froude*6.0*pi/0.34;
kx      = 2*pi/Lx;
dx      = Lx/nx;
x       = (0:nx-1)*dx;

%% FILENAME

fname = string.empty;
maxs = length(svec);
for s = 1:maxs
    fname(s) = sprintf('../%s/%s/%s/%s_s%d.h5', folder_name, data_folder, file_name, file_name, svec(s));
end

%% HDF5 FILE INFO

groupname = '/tasks';
info      = h5info(fname(1), groupname);
tasknames = {info.Datasets.Name};

%% DATA STORAGE PREALLOCATION

nvars   = length(tasknames);
gmodes  = struct('vars', repmat({zeros(nz,1)+1i*ones(nz,1)}, nvars, 1));

%% GET DATA FROM FILE

t = [];

for s = 1:maxs
   fprintf('Reading file %d of %d.\n',svec(s),svec(maxs)) 
   t = [t; h5read(fname(s),'/scales/sim_time')];
   z = h5read(fname(s),'/scales/z/1.0');
   for varindex = 1:nvars
       varname  = strcat(groupname,'/',tasknames{varindex});   
       tempvar  = h5read(fname(s),varname);
       tempvar2 = tempvar.r+1i*tempvar.i;
       if varindex == 5
           tempvar2 = repmat(tempvar2, nz, 1);
       end
       gmodes(varindex).vars = cat(2, gmodes(varindex).vars, tempvar2);
   end
end

%% DETERMINE TIMESERIES LENGTH

nf = length(t);

%% CLEAR OUT FIRST INDEX, RESHAPE, ADD x SHAPE FCN

modeindex = 1;
for varindex = 1:nvars
    gmodes(varindex).vars(:,1) = [];
    gmodes(varindex).vars = reshape(gmodes(varindex).vars, nz, 1, nf);
    gmodes(varindex).shapefcn = exp(1i*modes(modeindex)*kx*x);
    modeindex = modeindex + 1;
    if modeindex > length(modes)
        modeindex = 1;
    end
end

%% CALCULATE OUTPUT FIELDS

u    = zeros(nz,nx,nf);
w    = zeros(nz,nx,nf);
b    = zeros(nz,nx,nf);
vort = zeros(nz,nx,nf);
varstart = (0:3)*length(modes);

for i = 1:length(modes)
    b = b + pagemtimes(gmodes(varstart(1)+i).vars,gmodes(varstart(1)+i).shapefcn) ...
        + conj(pagemtimes(gmodes(varstart(1)+i).vars, gmodes(varstart(1)+i).shapefcn));
    u = u + pagemtimes(gmodes(varstart(2)+i).vars,gmodes(varstart(2)+i).shapefcn) ...
        + conj(pagemtimes(gmodes(varstart(2)+i).vars, gmodes(varstart(2)+i).shapefcn));
    w = w + pagemtimes(gmodes(varstart(3)+i).vars,gmodes(varstart(3)+i).shapefcn) ...
        + conj(pagemtimes(gmodes(varstart(3)+i).vars, gmodes(varstart(3)+i).shapefcn));
    vort = vort + pagemtimes(gmodes(varstart(4)+i).vars,gmodes(varstart(4)+i).shapefcn) ...
        + conj(pagemtimes(gmodes(varstart(4)+i).vars, gmodes(varstart(4)+i).shapefcn));
end

%% WRAP

if wrap == 1
    x    = cat(2, x, 2*x(end)-x(end-1));
    u    = cat(2, u, u(:, 1, :));
    w    = cat(2, w, w(:, 1, :));
    b    = cat(2, b, b(:, 1, :));
    vort = cat(2, vort, vort(:, 1, :));

    z    = cat(1, z, 2*z(end)-z(end-1));
    u    = cat(1, u, u(1, :, :));
    w    = cat(1, w, w(1, :, :));
    b    = cat(1, b, b(1, :, :));
    vort = cat(1, vort, vort(1, :, :));
end

end