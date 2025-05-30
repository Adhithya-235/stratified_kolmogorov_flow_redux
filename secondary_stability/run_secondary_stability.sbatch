#!/bin/bash
#
#SBATCH --job-name=secstab
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --output=slurm_logs/secstab_%A_%a.out
#SBATCH --error=slurm_logs/secstab_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as1751@usnh.edu

module purge
module load matlab

cd /mnt/lustre/space/as1751/Command_Center/Research/stratified_kolmogorov_flow_redux/secondary_stability

mkdir -p slurm_logs
mkdir -p solutions  # Ensure base folder exists safely before MATLAB runs

# Set constants
ALPHA=0
PR=1
FR=0.01
BETA=0
NUMEIG=30
TOL=1e-12
MAXIT=50000
VERBOSE=1

# Build a filtered list of ECS files that have not already been processed
ECSFILES=()
for file in import/ECS_real_field_for_Reb=*.mat; do
  reb=$(echo "$file" | sed -E 's/.*Reb=([0-9]+)\.mat/\1/')
  folder="solutions/Reb${reb}_alpha${ALPHA}"
  if [[ ! -d "$folder" ]]; then
    ECSFILES+=("$file")
  fi
done

# Gracefully exit if too few files
if [ $SLURM_PROCID -ge ${#ECSFILES[@]} ]; then
  echo "No ECS file assigned to proc $SLURM_PROCID. Skipping."
  exit 0
fi

ECSFILE=${ECSFILES[$SLURM_PROCID]}
REB=$(basename "$ECSFILE" | sed -E 's/[^0-9]*([0-9]+)\.mat/\1/')

# Call MATLAB with all values passed into the workspace
matlab -batch "clear; close all; clc; Reb=$REB; Pr=$PR; Fr=$FR; alpha=$ALPHA; beta=$BETA; numEig=$NUMEIG; residualtol=$TOL; maxit=$MAXIT; outputflag=$VERBOSE; test_ecs"

# Move SLURM logs into the corresponding output directory
SOLUTION_DIR=solutions/Reb${REB}_alpha${ALPHA}
mv slurm_logs/secstab_${SLURM_JOB_ID}_${SLURM_PROCID}.out "$SOLUTION_DIR/"
mv slurm_logs/secstab_${SLURM_JOB_ID}_${SLURM_PROCID}.err "$SOLUTION_DIR/"
