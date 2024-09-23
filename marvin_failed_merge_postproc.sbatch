#!/bin/bash
#SBATCH --job-name=sss_code_postproc
#SBATCH --output=/mnt/lustre/space/as1751/Command_Center/Research/stacked_stratified_shear/logger_postproc.out
#SBATCH --error=/mnt/lustre/space/as1751/Command_Center/Research/stacked_stratified_shear/logger_postproc.err
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
folder_name=2024-08-13_06-36-57
data_folder=results_mode2

# Make plotting subfolders
# mkdir -p "$folder_name/plots/timeseries"
# mkdir -p "$folder_name/plots/frames"
# mkdir -p "$folder_name/plots/movies"
# mkdir -p "$folder_name/plots/tke_hovmollers"

# Merge processes
mpiexec -n 24 python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/field_snapshots
mpiexec -n 24 python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/energy_timeseries
mpiexec -n 24 python3 -m dedalus merge_procs "$folder_name"/"$data_folder"/checkpointing_data

# Move generated data and logging files to the created folder
mv logger_postproc.* "$folder_name/"
cp marvin_failed_merge_postproc.sbatch "$folder_name/"

# Create analysis plots after navigating to post processing code folder
cd $HOME/Command_Center/Research/stacked_stratified_shear/post_processing/
echo "Plotting KE timeseries."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; maxs=1; plot_te"
echo "Plotting probe frames."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=10; svec=[1:20]; wrap=0; plot_fields"
echo "Plotting perturbation energies."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:20]; wrap=1; unwrap=0; Ri=0.2; plot_pte"
echo "Plotting timeseries of buoyancy power spectrum."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:20]; wrap=1; unwrap=1;  plot_spectrumtime"
echo "Plotting interfacial buoyancy hovmoller plot."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:20]; wrap=1; unwrap=0;  plot_spacetime"
echo "Performing TKE budget analysis."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:20]; wrap=1; unwrap=0; Ri=0.2; Re=10000; plot_tke_budget"
echo "Animating buoyancy slices."
matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$folder_name'; data_folder='$data_folder'; file_name='field_snapshots'; stride=1; svec=[1:20]; wrap=0; animate_fields"

# Display a message indicating completion
echo "Simulation complete. Post processing complete. Results saved in $folder_name."