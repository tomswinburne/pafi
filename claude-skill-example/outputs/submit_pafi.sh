#!/bin/bash
#SBATCH --job-name=pafi-w-vac
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --account=mat-w-pafi

# --- cluster environment (preserved verbatim by the skill) -----------------
module purge
module load gcc/12.2.0 openmpi/4.1.5 lammps/2024.06.27
source ~/venvs/pafi/bin/activate

cd $SLURM_SUBMIT_DIR

# --- launch line (the skill replaces <LAUNCH>) -----------------------------
mpirun -n $SLURM_NTASKS python run_pafi.py
