#!/bin/bash
#SBATCH --job-name=fe_NO2
#SBATCH --output=qe_job.out
#SBATCH --error=qe_job.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=04:00:00
#SBATCH --partition=standard

module load quantum-espresso

mpirun -np 16 pw.x -in qespresso_input.in > qe_output.out
