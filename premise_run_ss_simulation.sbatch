#!/bin/bash
#SBATCH --job-name=sskf_code
#SBATCH --output=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stratified_kolmogorov_flow_redux/logger_branch1.out
#SBATCH --error=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stratified_kolmogorov_flow_redux/logger_branch1.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
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

# Create a folder with the current date and time
folder_name=$(date +'%Y-%m-%d_%H-%M-%S')
mkdir "$folder_name"
data_folder=results_branch1
mkdir "$data_folder"

# Make plotting subfolders
mkdir -p "$folder_name/plots/timeseries"
mkdir -p "$folder_name/plots/frames"
mkdir -p "$folder_name/plots/movies"

# Run simulation, merge processes
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 stratified_kolmogorov_flow.py --Rb=50.0 --Pr=1.0 --Fr=0.02 --Nx=1024 --Nz=1024 --Tend=15.0 
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$data_folder"/field_snapshots
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$data_folder"/energy_timeseries
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$data_folder"/checkpointing_data

# Move generated data and logging files to the created folder
mv "$data_folder/" logger_branch1.* "$folder_name/"
cp stratified_kolmogorov_flow.py premise_run_ss_simulation.sbatch "$folder_name/"

# Create analysis plots after navigating to post processing code folder
cd $HOME/Research_Active_Projects/stratified_kolmogorov_flow_redux/post_processing/
echo "Plotting KE timeseries."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; maxs=1; plot_te"
echo "Plotting probe frames."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=100; svec=[1:15]; wrap=0; Fr = 0.02; plot_fields"
echo "Plotting perturbation energies."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:15]; wrap=1; unwrap=0; Fr=0.02; plot_pte"
echo "Plotting spacetime vorticity."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:15]; wrap=1; unwrap=0; Fr=0.02; plot_spacetime"
echo "Animating vorticity slices."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:15]; wrap=0; Fr = 0.02; animate_fields"

# Display a message indicating completion
echo "Simulation complete. Post processing complete. Results saved in $folder_name."