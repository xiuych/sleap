#!/bin/bash

#SBATCH --job-name=sleap
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-22500
#SBATCH --output=/proj/ibrahimlab/leap/sims_pwe3/slurm/_sleap/slurmLogFiles_sleap_%a.out
#SBATCH --error=/proj/ibrahimlab/leap/sims_pwe3/slurm/_sleap/sleap_%a.err
#SBATCH --constraint=rhel8

## add R module
# module add gcc/6.3.0
module add r/4.1.0

R CMD BATCH --no-restore /proj/ibrahimlab/leap/sims_pwe3/02_run_sims.R /proj/ibrahimlab/leap/sims_pwe3/slurm/_sleap/sleap_$SLURM_ARRAY_TASK_ID.Rout