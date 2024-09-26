#!/bin/bash
#SBATCH --job-name=sskf_code
#SBATCH --output=/mnt/lustre/space/as1751/Command_Center/Research/stratified_kolmogorov_flow_redux/logger_postproc.out
#SBATCH --error=/mnt/lustre/space/as1751/Command_Center/Research/stratified_kolmogorov_flow_redux/logger_postproc.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --time=9000

# Load module and set other things
export OMP_NUM_THREADS=1 
module load Anaconda3/2022.05
module load MATLAB

# Activate dedalus
source /mnt/lustre/software/Anaconda3/Anaconda3-2022.05/etc/profile.d/conda.sh
conda activate dedalus-20221006

# Point to solution folder system
cd $HOME/Command_Center/Research/stacked_stratified_shear/
folder_name=2024-09-26_10-38-56
data_folder=results_branch4

# Merge processes
mpiexec -n 24 python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/field_snapshots
mpiexec -n 24 python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/energy_timeseries
mpiexec -n 24 python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/checkpointing_data

# Move generated data and logging files to the created folder
mv logger_postproc.* "$folder_name/"
cp marvin_failed_merge_postproc.sbatch "$folder_name/"

# Create analysis plots after navigating to post processing code folder
cd $HOME/Command_Center/Research/stratified_kolmogorov_flow_redux/post_processing/
echo "Plotting KE timeseries."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; maxs=1; plot_te"
echo "Plotting probe frames."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=100; svec=[1:100]; wrap=0; Fr = 0.02; plot_fields"
echo "Plotting perturbation energies."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:100]; wrap=1; unwrap=0; Fr=0.02; plot_pte"
echo "Animating vorticity slices."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:100]; wrap=0; Fr = 0.02; animate_fields"

# Display a message indicating completion
echo "Simulation complete. Post processing complete. Results saved in $folder_name."