#!/bin/bash
#SBATCH --job-name=sskf_code
#SBATCH --output=/mnt/lustre/space/as1751/Command_Center/Research/stratified_kolmogorov_flow_redux/logger_branch3.out
#SBATCH --error=/mnt/lustre/space/as1751/Command_Center/Research/stratified_kolmogorov_flow_redux/logger_branch3.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --time=9000

# Load module and set other things
export OMP_NUM_THREADS=1 
module load Anaconda3/2022.05
module load MATLAB

# Activate dedalus
source /mnt/lustre/software/Anaconda3/Anaconda3-2022.05/etc/profile.d/conda.sh
conda activate dedalus-20221006

# Create a folder with the current date and time
cd $HOME/Command_Center/Research/stratified_kolmogorov_flow_redux/
folder_name=$(date +'%Y-%m-%d_%H-%M-%S')
mkdir "$folder_name"
data_folder=results_branch3
mkdir "$data_folder"

# Make plotting subfolders
mkdir -p "$folder_name/plots/timeseries"
mkdir -p "$folder_name/plots/frames"
mkdir -p "$folder_name/plots/movies"

# Run simulation, merge processes
mpiexec -n 16 python3 ~/Command_Center/Research/stratified_kolmogorov_flow_redux/stratified_kolmogorov_flow.py --Rb=50.0 --Pr=1.0 --Fr=0.02 --Nx=512 --Nz=512 --Tend=50.0 
mpiexec -n 16 python3 -m dedalus merge_procs "$data_folder"/field_snapshots
mpiexec -n 16 python3 -m dedalus merge_procs "$data_folder"/energy_timeseries
mpiexec -n 16 python3 -m dedalus merge_procs "$data_folder"/checkpointing_data

# Move generated data and logging files to the created folder
mv "$data_folder/" logger_branch3.* "$folder_name/"
cp stratified_kolmogorov_flow.py marvin_run_ss_simulation.sbatch "$folder_name/"

# Create analysis plots after navigating to post processing code folder
cd $HOME/Command_Center/Research/stratified_kolmogorov_flow_redux/post_processing/
echo "Plotting KE timeseries."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; maxs=1; plot_te"
echo "Plotting probe frames."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=100; svec=[1:50]; wrap=0; Fr = 0.02; plot_fields"
echo "Plotting perturbation energies."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:50]; wrap=1; unwrap=0; Fr=0.02; plot_pte"
echo "Animating vorticity slices."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:50]; wrap=0; Fr = 0.02; animate_fields"

# Display a message indicating completion
echo "Simulation complete. Post processing complete. Results saved in $folder_name."