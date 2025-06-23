#!/bin/bash

#SBATCH --job-name=gendata
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000m
#SBATCH --array=1-45
#SBATCH --output=/proj/ibrahimlab/leap/sims_pwe3/slurm/slurmLogFiles_gendata_%a.out
#SBATCH --error=/proj/ibrahimlab/leap/sims_pwe3/slurm/gendata_%a.err
#SBATCH --constraint=rhel8

## add R module
# module add gcc/6.3.0
module add r/4.1.0

R CMD BATCH --no-restore /proj/ibrahimlab/leap/sims_pwe3/01_generate_data.R /proj/ibrahimlab/leap/sims_pwe3/slurm/gendata_$SLURM_ARRAY_TASK_ID.Rout