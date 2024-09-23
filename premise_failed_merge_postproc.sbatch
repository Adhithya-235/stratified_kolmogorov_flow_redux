#!/bin/bash
#SBATCH --job-name=sss_code_postproc
#SBATCH --output=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stacked_stratified_shear/logger_postproc.out
#SBATCH --error=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stacked_stratified_shear/logger_postproc.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=4
#SBATCH --time=9000

# Setup 

export OMP_NUM_THREADS=1 
export NUMEXPR_MAX_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

# Navigate to code directory, load module 
cd $HOME/Research_Active_Projects/stacked_stratified_shear/
module load Anaconda3/2019.10
module load MATLAB/R2022b_Update_3

# Point to solution folder system
folder_name=2024-08-05_23-36-49
data_folder=results_twomode

# Merge processes
mpiexec -n 96 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/field_snapshots
mpiexec -n 96 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/energy_timeseries
mpiexec -n 96 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/checkpointing_data

# Move generated data and logging files to the created folder
mv logger_postproc.* "$folder_name/"
cp failed_merge_postproc.sbatch "$folder_name/"

# Create analysis plots after navigating to post processing code folder
cd $HOME/Research_Active_Projects/stacked_stratified_shear/post_processing/
echo "Plotting KE timeseries."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; maxs=1; plot_te"
echo "Plotting probe frames."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=10; svec=[1:19]; wrap=0; plot_fields"
echo "Plotting perturbation energies."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:19]; wrap=1; unwrap=0; Ri=0.2; plot_pte"
echo "Plotting timeseries of buoyancy power spectrum."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:19]; wrap=1; unwrap=1;  plot_spectrumtime"
echo "Plotting interfacial buoyancy hovmoller plot."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:19]; wrap=1; unwrap=0;  plot_spacetime"
echo "Performing TKE budget analysis."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:20]; wrap=1; unwrap=0; Ri=0.2; Re=10000; plot_tke_budget"
echo "Animating buoyancy slices."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:19]; wrap=0; animate_fields"

# Display a message indicating completion
echo "Post processing complete. Results saved in $folder_name."