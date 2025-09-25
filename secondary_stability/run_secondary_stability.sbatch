#!/bin/bash
#SBATCH --job-name=secstab
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00
#SBATCH --chdir=/mnt/lustre/space/as1751/Command_Center/Research/stratified_kolmogorov_flow_redux/secondary_stability
#SBATCH --mail-user=as1751@usnh.edu
#SBATCH --mail-type=ALL

set -euo pipefail
module load MATLAB

# 1) Prepare log + solution roots

mkdir -p slurm_logs solutions

# 2) Gather all .mat ECS files into a Bash array

ECSFILES=(import/ECS_real_field_for_Reb=*.mat)

# 3) Optional sanity check: fail early if file count ≠ cores requested

if [[ ${#ECSFILES[@]} -ne $SLURM_NTASKS ]]; then
  echo "Found ${#ECSFILES[@]} files but SLURM_NTASKS=$SLURM_NTASKS" >&2
  exit 1
fi

# 4) Launch NTASKS parallel tasks, one per core

srun \
 --ntasks=$SLURM_NTASKS \
 --cpus-per-task=1 \
 --output=slurm_logs/secstab_%j_%t.out \
 --error=slurm_logs/secstab_%j_%t.err \
 bash -c '

  # a) Pick the file for this task

  idx=$(( SLURM_PROCID + 1 ))
  FILE="${!idx}"
  echo "[$SLURM_PROCID] Processing file: $FILE" >&2
  REB_VAL=$(basename "$FILE" | sed -E "s/[^0-9]*([0-9]+)\.mat/\1/")
  echo "[$SLURM_PROCID] Parsed Reb = $REB_VAL" >&2
  ALPHA_VAL=0.5

  # b) Format names exactly as the MATLAB script does

  REB=$(printf "%.2f" "$REB_VAL")
  ALPHA=$(printf "%.2f" "$ALPHA_VAL")
  solution_dir="solutions/Reb${REB}_alpha${ALPHA}"

  # c) Skip if MATLAB has already created that folder

  if [[ -d "$solution_dir" ]]; then

    echo "[$SLURM_PROCID] Skipping Reb=${REB_VAL} (solution exists)" >&2

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

' bash "${ECSFILES[@]}"
