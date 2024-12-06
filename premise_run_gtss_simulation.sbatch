#!/bin/bash
#SBATCH --job-name=sskf_code_galerkin
#SBATCH --output=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stratified_kolmogorov_flow_redux/logger.out
#SBATCH --error=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stratified_kolmogorov_flow_redux/logger.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
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
data_folder=results_galerkin
mkdir "$data_folder"

# Make plotting subfolders
mkdir -p "$folder_name/plots/timeseries"
mkdir -p "$folder_name/plots/frames"
mkdir -p "$folder_name/plots/movies"

# Run simulation, merge processes
mpiexec -n 1 $HOME/.conda/envs/dedalus/bin/python3 galerkin_kolmogorov_flow.py --Rb=50.0 --Pr=1.0 --Fr=0.02 --Nz=512 --Tend=15.0 
mpiexec -n 1 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs "$data_folder"/field_snapshots

# Move generated data and logging files to the created folder
mv "$data_folder/" logger.* "$folder_name/"
cp galerkin_kolmogorov_flow.py premise_run_gtss_simulation.sbatch "$folder_name/"

# Create analysis plots after navigating to post processing code folder
cd $HOME/Research_Active_Projects/stratified_kolmogorov_flow_redux/post_processing_galerkin/
echo "Plotting perturbation energies."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; nx=512; nz=512; svec=[1:15]; wrap=1; unwrap=0; Fr=0.02; plot_pte"
echo "Animating vorticity slices."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; nx=512; nz=512; svec=[1:15]; wrap=0; Fr = 0.02; animate_fields"


# Display a message indicating completion
echo "Simulation complete. Post processing complete. Results saved in $folder_name."