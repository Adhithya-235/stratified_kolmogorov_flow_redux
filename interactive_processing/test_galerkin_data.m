clear
close all
clc

%% ADD PATH

addpath('../utility_belt'); 

%% FOLDER DETAILS

folder_name='2024-12-09_09-57-04'; 
data_folder='results_galerkin'; 
file_name='field_snapshots'; 

%% RUNTIME PARAMETERS

nx=1024; 
nz=1024; 
svec=1:4; 
modes=[0, 19, 23, 4]; 
wrap=0; 
unwrap=0; 
Fr=0.02; 

%% READ DATA

[x, z, t, u, w, b, ~, nf] = get_galerkin_data(folder_name, data_folder, file_name, svec, nx, nz, modes, wrap);