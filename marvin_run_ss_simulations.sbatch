#!/bin/bash
#SBATCH --job-name=skf_sweep
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1       # Request all tasks on a single node
#SBATCH --ntasks=80     # Request 10 sims * 8 CPUs/sim = 80 CPUs total
#SBATCH --cpus-per-task=1 # Each task gets 1 CPU
#SBATCH --time=9000
#
# Main job log (slurm-JOB_ID.out) will show the overall status.
# Each individual simulation's output will be in:
# symmetrybreak/Reb.../run.log

#====================================================================
# === 1. DEFINE THE SIMULATION FUNCTION ===
# This function runs ONE simulation
# Usage: run_simulation <RB_VALUE> <EFNID_VALUE>
#====================================================================
run_simulation() {
    # --- Assign arguments to local variables ---
    # Using 'local' prevents variables from interfering
    # with other parallel runs.
    local RB=$1
    local EFNID=$2

    # --- Set other parameters (constant for all runs) ---
    local PR=1.00
    local FR=0.01
    local ALPHA=0.00
    local NX=128
    local NZ=128
    local TEND=50.0
    local N_PROCS_PER_SIM=8 # CPUs for this one simulation

    # --- Generate paths for this specific run ---
    local RUN_ID=$(printf "Reb%.2f_alpha%.2f_idx%d" $RB $ALPHA $EFNID)
    local RUN_DIR="symmetrybreak/$RUN_ID"
    local LOG_FILE="$RUN_DIR/run.log"
    
    # --- Create directories and archive scripts ---
    mkdir -p "$RUN_DIR/plots/timeseries"
    mkdir -p "$RUN_DIR/plots/frames"
    mkdir -p "$RUN_DIR/plots/movies"
    
    # $0 is the name of this sbatch file
    cp stratified_kolmogorov_flow.py "$0" "$RUN_DIR/"

    echo "--- Starting simulation: $RUN_ID ---"
    echo "Running Dedalus... Output will be in $LOG_FILE"

    # --- 4. RUN SIMULATION ---
    # Use 'srun' to launch this 8-task job.
    # --exclusive ensures it gets its 8 CPUs without sharing.
    # Redirect all stdout/stderr to this run's specific log file.
    mpiexec -n $N_PROCS_PER_SIM \
    python3 stratified_kolmogorov_flow.py \
        --Rb=$RB \
        --Pr=$PR \
        --Fr=$FR \
        --Nx=$NX \
        --Nz=$NZ \
        --Tend=$TEND \
        --Efnid=$EFNID \
        --alpha=$ALPHA > "$LOG_FILE" 2>&1

    echo "Merging simulation output..." >> "$LOG_FILE"

    # --- 5. MERGE DATA ---
    # Use 'srun' again for the merging steps.
    # Append output to the same log file.
    mpiexec -n $N_PROCS_PER_SIM \
    python3 -m dedalus merge_procs "$RUN_DIR"/field_snapshots --cleanup >> "$LOG_FILE" 2>&1
    
    mpiexec -n $N_PROCS_PER_SIM \
    python3 -m dedalus merge_procs "$RUN_DIR"/energy_timeseries --cleanup >> "$LOG_FILE" 2>&1
    
    mpiexec -n $N_PROCS_PER_SIM \
    python3 -m dedalus merge_procs "$RUN_DIR"/checkpointing_data --cleanup >> "$LOG_FILE" 2>&1

    # --- 6. RUN POST-PROCESSING ---
    echo "Running MATLAB post-processing..." >> "$LOG_FILE"
    
    # Define paths for MATLAB
    local SCRIPT_DIR="$HOME/Command_Center/Research/stratified_kolmogorov_flow_redux"
    local POST_DIR="$SCRIPT_DIR/post_processing"
    local FOLDER_NAME_MATLAB="symmetrybreak"
    local DATA_FOLDER_MATLAB="$RUN_ID"

    # Run matlab in a subshell (...) to prevent 'cd' from
    # affecting other parallel runs. Append output to log.
    
    # TASK 1: Plot Energy Timeseries
    echo "Plotting KE timeseries..." >> "$LOG_FILE"
    (cd "$POST_DIR" && \
    matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$FOLDER_NAME_MATLAB'; data_folder='$DATA_FOLDER_MATLAB'; maxs=1; plot_te") >> "$LOG_FILE" 2>&1

    # TASK 2: Plot Probe Frames
    echo "Plotting probe frames..." >> "$LOG_FILE"
    (cd "$POST_DIR" && \
    matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$FOLDER_NAME_MATLAB'; data_folder='$DATA_FOLDER_MATLAB'; file_name='field_snapshots'; stride=100; svec=[1:$TEND]; wrap=0; Fr = $FR; plot_fields") >> "$LOG_FILE" 2>&1
    
    # TASK 3: Calculate Perturbation Energies
    echo "Calculate perturbation energies ..." >> "$LOG_FILE"
    (cd "$POST_DIR" && \
    matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$FOLDER_NAME_MATLAB'; data_folder='$DATA_FOLDER_MATLAB'; file_name='field_snapshots'; stride=1; svec=[1:$TEND]; Nx = $NX; Nz = $NZ; alpha = $ALPHA; plot_secstab_pte") >> "$LOG_FILE" 2>&1

    # TASK 4: Animate Vorticity Slices
    echo "Animating vorticity slices..." >> "$LOG_FILE"
    (cd "$POST_DIR" && \
    matlab -batch "clear; close all; clc; addpath('../utility_belt'); folder_name='$FOLDER_NAME_MATLAB'; data_folder='$DATA_FOLDER_MATLAB'; file_name='field_snapshots'; stride=1; svec=[1:$TEND]; wrap=0; Fr = $FR; animate_fields") >> "$LOG_FILE" 2>&1
    
    echo "--- Finished simulation: $RUN_ID ---" >> "$LOG_FILE"
}

#====================================================================
# === 2. SET UP MAIN ENVIRONMENT (Done once) ===
#====================================================================
echo "Job $SLURM_JOB_ID starting."
echo "Requesting $SLURM_NTASKS CPUs on 1 node."

# Set OMP threads
export OMP_NUM_THREADS=1 

# Load modules
module load Anaconda3/2022.05
module load MATLAB

# Activate dedalus environment
source /mnt/lustre/software/Anaconda3/Anaconda3-2022.05/etc/profile.d/conda.sh
conda activate dedalus-20221006

# Navigate to the code directory
cd $HOME/Command_Center/Research/stratified_kolmogorov_flow_redux/

# Make the function available to srun
export -f run_simulation

#====================================================================
# === 3. LAUNCH ALL SIMULATIONS IN PARALLEL ===
#====================================================================
echo "Launching 10 simulations in parallel..."

# Rbs 4-14 with Efnid=1
for RB_VAL in 4 6 8 10 12 14; do
    run_simulation $RB_VAL 1 &
done

# Rb 16 with Efnid=1, 2
run_simulation 16 1 &
run_simulation 16 2 &

# Rb 18 with Efnid=1, 3
run_simulation 18 1 &
run_simulation 18 3 &

#====================================================================
# === 4. WAIT ===
#====================================================================
# 'wait' tells the script to pause here until all
# background jobs (&) have finished.
echo "All 10 simulations launched. Waiting for completion..."
wait

echo "All simulations finished. Job $SLURM_JOB_ID complete."