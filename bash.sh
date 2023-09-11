#!/bin/bash
#SBATCH --job-name=a.out
#SBATCH --partition=mono
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
##SBATCH --cpus-per-task=1
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --time=0-02:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
. /etc/profile

./a.out
