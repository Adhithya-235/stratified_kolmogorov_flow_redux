function [t, nf, ke, pe, te] = get_timeseries_data(data_folder,folder_name, maxs)

% This function reads HDF5 data produced by dedalus and extracts temporal
% grid as well as timeseries of the total and cross-wind kinetic energies. 
% Specify date-based folder name as well as the maximum number of data series produced. 
% Also, maxs is the maximum series number.
    
%% FILENAME
    
fname = string.empty;
    
for s = 1:maxs
    fname(s) = sprintf('../%s/%s/energy_timeseries/energy_timeseries_s%d.h5', folder_name, data_folder, s);
end

%% READ DATA
    
t = [];
ke = [];
pe = [];
te = [];

for s = 1:maxs
    t = [t; h5read(fname(s),'/scales/sim_time')];
    tempke = h5read(fname(s),'/tasks/KE');
    temppe = h5read(fname(s),'/tasks/PE');
    tempte = h5read(fname(s),'/tasks/TE');
    ke = [ke; tempke(:)];
    pe = [pe; temppe(:)];
    te = [te; tempte(:)];
end

nf = length(t);

end