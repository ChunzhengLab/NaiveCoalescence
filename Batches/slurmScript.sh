#!/bin/bash
#SBATCH --partition=LongTermPool
#SBATCH --job-name=Coalnt
#SBATCH --array=0-15
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=output.log

srun ./commands.sh $SLURM_ARRAY_TASK_ID
