#!/bin/bash
#SBATCH --job-name=<JOBNAME>
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=<NTASKS>
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --account=mat-w-pafi

# --- cluster environment (preserved verbatim by the skill) -----------------
module purge
module load gcc/12.2.0 openmpi/4.1.5 lammps/2024.06.27
source ~/venvs/pafi-ase/bin/activate

cd $SLURM_SUBMIT_DIR

# --- launch line (the skill replaces <LAUNCH>) -----------------------------
<LAUNCH>
