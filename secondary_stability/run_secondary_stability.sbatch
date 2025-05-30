#!/bin/bash
#SBATCH --job-name=secstab
#SBATCH --array=0-19                # one task per ECS file
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --chdir=/mnt/lustre/space/as1751/Command_Center/Research/stratified_kolmogorov_flow_redux/secondary_stability
#SBATCH --output=slurm_logs/secstab_%A_%a.out
#SBATCH --error=slurm_logs/secstab_%A_%a.err
#SBATCH --mail-user=as1751@usnh.edu
#SBATCH --mail-type=FAIL,END

set -euo pipefail
module purge
module load MATLAB

# make sure the log and result directories exist *before* anything else
mkdir -p slurm_logs solutions

# Select the ECS file that corresponds to this array index
ECSFILES=(import/ECS_real_field_for_Reb=*.mat)
if [[ ${SLURM_ARRAY_TASK_ID} -ge ${#ECSFILES[@]} ]]; then
    echo "Nothing to do for index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi
ECSFILE=${ECSFILES[$SLURM_ARRAY_TASK_ID]}
REB=$(basename "$ECSFILE" | sed -E 's/[^0-9]*([0-9]+)\.mat/\1/')

# Parameters
PR=1 ; FR=0.01 ; ALPHA=0 ; BETA=0
NUMEIG=30 ; TOL=1e-12 ; MAXIT=50000 ; VERBOSE=1

matlab -batch "clear;close all;clc; Reb=$REB; Pr=$PR; Fr=$FR; \
               alpha=$ALPHA; beta=$BETA; numEig=$NUMEIG; \
               residualtol=$TOL; maxit=$MAXIT; outputflag=$VERBOSE; \
               ecs_stability"

# organise log files
mv slurm_logs/secstab_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out  \
   solutions/Reb${REB}_alpha${ALPHA}/
mv slurm_logs/secstab_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err  \
   solutions/Reb${REB}_alpha${ALPHA}/
