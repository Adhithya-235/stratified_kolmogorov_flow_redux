#!/bin/bash
#SBATCH --job-name=sskf_code_postproc
#SBATCH --output=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stratified_kolmogorov_flow_redux/logger_postproc.out
#SBATCH --error=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stratified_kolmogorov_flow_redux/logger_postproc.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=9000

# Setup 

export OMP_NUM_THREADS=1 
export NUMEXPR_MAX_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

# Navigate to code directory, load module 
cd $HOME/Research_Active_Projects/stratified_kolmogorov_flow_redux/
module load Anaconda3/2019.10
module load MATLAB/R2022b_Update_3

# Point to solution folder system
folder_name=2024-09-23_09-09-44
data_folder=results_branch4

# Merge processes
mpiexec -n 24 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/field_snapshots
mpiexec -n 24 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/energy_timeseries
mpiexec -n 24 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/checkpointing_data

# Move generated data and logging files to the created folder
mv logger_postproc.* "$folder_name/"
cp premise_failed_merge_postproc.sbatch "$folder_name/"

# Create analysis plots after navigating to post processing code folder
cd $HOME/Research_Active_Projects/stratified_kolmogorov_flow_redux/post_processing/
echo "Plotting KE timeseries."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; maxs=1; plot_te"
echo "Plotting probe frames."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=100; svec=[1:100]; wrap=0; Fr = 0.02; plot_fields"
echo "Plotting perturbation energies."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:100]; wrap=1; unwrap=0; Fr=0.02; plot_pte"
echo "Animating vorticity slices."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:100]; wrap=0; Fr = 0.02; animate_fields"

# Display a message indicating completion
echo "Post processing complete. Results saved in $folder_name."