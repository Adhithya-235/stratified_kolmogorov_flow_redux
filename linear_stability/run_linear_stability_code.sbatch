#!/bin/bash
#SBATCH --job-name=viscstab_code
#SBATCH --output=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stacked_stratified_shear/viscous_stability_boundaries/logger.out
#SBATCH --error=/mnt/lz01/chini/as1751-new/Research_Active_Projects/stacked_stratified_shear/viscous_stability_boundaries/logger.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=9000

# Navigate to code directory, load module 
cd $HOME/Research_Active_Projects/stacked_stratified_shear/viscous_stability_boundaries/
module load MATLAB/R2022b_Update_3

# Create a folder with the current date and time, copy current version of code
folder_name=$(date +'%Y-%m-%d_%H-%M-%S')
mkdir "$folder_name"
cp plot_stability_bndry.m "$folder_name/"

# Run simulation
echo "Plotting stability boundaries."
matlab -batch "plot_stability_bndry"

# Move generated data and logging files to the created folder
mv *.mat *.fig *.png logger.* "$folder_name/"

# Display a message indicating completion
echo "Simulation complete. Post processing complete. Results saved in $folder_name."