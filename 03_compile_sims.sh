#!/bin/bash

#SBATCH --job-name=compileSims
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000m
#SBATCH --array=1-45
#SBATCH --output=/proj/ibrahimlab/leap/sims_pwe3/slurm/slurmLogFiles_compileSims_%a.out
#SBATCH --error=/proj/ibrahimlab/leap/sims_pwe3/slurm/compileSims_%a.err
#SBATCH --constraint=rhel8

## add R module
# module add gcc/6.3.0
module add r/4.1.0

R CMD BATCH --no-restore /proj/ibrahimlab/leap/sims_pwe3/03_compile_sims.R /proj/ibrahimlab/leap/sims_pwe3/slurm/compileSims_$SLURM_ARRAY_TASK_ID.Rout