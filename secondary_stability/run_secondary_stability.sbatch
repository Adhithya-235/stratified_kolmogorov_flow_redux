#!/bin/bash
#SBATCH --job-name=secstab_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=112            # <-- MODIFIED: 4 * 28 = 112 tasks
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00
#SBATCH --chdir=/mnt/lustre/space/as1751/Command_Center/Research/stratified_kolmogorov_flow_redux/secondary_stability
#SBATCH --mail-user=as1751@usnh.edu
#SBATCH --mail-type=ALL

set -euo pipefail
module load MATLAB

# 1) Prepare log + solution roots
mkdir -p slurm_logs solutions

# 2) Define all ALPHA values
ALPHA_VALUES=(0.50 0.3333333333333333 0.25 0.20) # Example 4 values of alpha

# 3) Gather all .mat ECS files into a Bash array
ECSFILES=(import/ECS_real_field_for_Reb=*.mat)
NUM_ECS_FILES=${#ECSFILES[@]}
NUM_ALPHA_VALUES=${#ALPHA_VALUES[@]}
EXPECTED_TASKS=$(( NUM_ECS_FILES * NUM_ALPHA_VALUES ))

# 4) Optional sanity check: fail early if file count * alpha count ≠ cores requested
if [[ $EXPECTED_TASKS -ne $SLURM_NTASKS ]]; then
  echo "Found ${NUM_ECS_FILES} files and ${NUM_ALPHA_VALUES} alpha values. Expected $EXPECTED_TASKS tasks but SLURM_NTASKS=$SLURM_NTASKS" >&2
  exit 1
fi

# 5) Generate a MASTER LIST of (Reb_val, Alpha_val) pairs
declare -a MASTER_TASK_LIST

for ALPHA_VAL_ITEM in "${ALPHA_VALUES[@]}"; do
  for ECS_FILE_ITEM in "${ECSFILES[@]}"; do
    # Append the Reb/Alpha pair as a single string to the master list
    REB_VAL=$(basename "$ECS_FILE_ITEM" | sed -E "s/[^0-9]*([0-9]+)\.mat/\1/")
    MASTER_TASK_LIST+=("${REB_VAL},${ALPHA_VAL_ITEM}")
  done
done

# 5a) Logging: Print a header and all generated task combinations
echo ""
echo "========================================================"
echo "LOGGING MASTER TASK LIST (Total Tasks: ${#MASTER_TASK_LIST[@]})"
echo "========================================================"

# Use printf to create a simple, aligned table header
printf "%-10s %-10s %-15s\n" "INDEX" "REB_VAL" "ALPHA_VAL"
printf "%-10s %-10s %-15s\n" "-----" "-------" "---------"

# 5b) Loop through the generated list and print each pair
for i in "${!MASTER_TASK_LIST[@]}"; do
  # The task pair is a string like "Reb,Alpha"
  TASK_PAIR="${MASTER_TASK_LIST[$i]}"
  
  # Split the "Reb,Alpha" string into two separate variables
  IFS=',' read -r REB_VAL ALPHA_VAL <<< "$TASK_PAIR"

  # Print the index (i), Reb value, and Alpha value
  printf "%-10s %-10s %-15s\n" "$i" "$REB_VAL" "$ALPHA_VAL"
done

echo "========================================================"
echo ""

# 6) Launch NTASKS parallel tasks, one per core, using the master list

srun \
 --ntasks=$SLURM_NTASKS \
 --cpus-per-task=1 \
 --output=slurm_logs/secstab_%j_%t.out \
 --error=slurm_logs/secstab_%j_%t.err \
 bash -c '

  # a) Pick the Reb and Alpha for this task from the MASTER_TASK_LIST
  # The MASTER_TASK_LIST is passed as trailing arguments to the second bash

  idx=$(( SLURM_PROCID + 1 ))
  TASK_PAIR="${!idx}"
  
  # Split the "Reb,Alpha" string into two separate variables
  IFS=',' read -r REB_VAL ALPHA_VAL <<< "$TASK_PAIR"

  echo "[$SLURM_PROCID] Processing Reb=${REB_VAL} with Alpha=${ALPHA_VAL}" >&2

  # b) Format names exactly as the MATLAB script does

  REB=$(printf "%.2f" "$REB_VAL")
  ALPHA=$(printf "%.2f" "$ALPHA_VAL")
  solution_dir="solutions/Reb${REB}_alpha${ALPHA}"

  # c) Skip if MATLAB has already created that folder

  if [[ -d "$solution_dir" ]]; then

    echo "[$SLURM_PROCID] Skipping Reb=${REB_VAL}, Alpha=${ALPHA_VAL} (solution exists)" >&2

  else

    # d) Run MATLAB; ecs_stability creates solution_dir itself

    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    matlab -singleCompThread -batch "clear; close all; clc; Reb=${REB_VAL}; Pr=1; Fr=0.01; alpha=${ALPHA_VAL}; beta=0; numEig=30; residualtol=1e-12; maxit=10000; outputflag=1; ecs_stability;"
  fi

  # e) Move this task’s Slurm logs into its solution folder

  mkdir -p "$solution_dir"
  mv slurm_logs/secstab_${SLURM_JOB_ID}_${SLURM_PROCID}.out "$solution_dir"/
  mv slurm_logs/secstab_${SLURM_JOB_ID}_${SLURM_PROCID}.err "$solution_dir"/

' bash "${MASTER_TASK_LIST[@]}"