"""
Example MPI usage of PAFIMPI with a LAMMPSlib EAM calculator.

Run with:
    mpirun -n 4 python mpi-ase-test.py

Each rank owns an independent LAMMPSlib/LAMMPS instance.  

All ranks sample the same (r, T) simultaneously.

Default assumes one MPI process/worker, required by LAMMPSLib in ASE

- inter_comm does the averaging / gathering (defaults to WORLD)
- worker_comm handles the local calculator (defaults to SELF)

In all cases, worker_comm MUST match the comm passed to the calculator

"""

import os
import numpy as np
from mpi4py import MPI
from ase.io import read
from ase.calculators.lammpslib import LAMMPSlib
from pafi import PAFI, PAFIMPI

# we will add to this file if it exists
json_dump_file = f"pafi_data_all_T_all_r.json"

# Default: 1 process per worker
worker_comm = MPI.COMM_SELF
inter_comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()

# Multi-CPU example (Not possible for LAMMPSLib):
# w_cpu = 2
# worker_comm = world.Split(color=rank//w_cpu, key=rank)
# inter_comm = world.Split(color=0 if rank%w_cpu==0 else MPI.UNDEFINED, key=rank//w_cpu)

# ---------------------------------------------------------------------------
# Load NEB images (all ranks read the same trajectory)
# ---------------------------------------------------------------------------
images = read("neb_paths/w-vac-path.traj", index=":")

# ---------------------------------------------------------------------------
# Each rank creates its own LAMMPSlib instance (one LAMMPS process per rank).
# keep_alive=True avoids restarting LAMMPS between force calls.
# ---------------------------------------------------------------------------
eam_path = "neb_paths/W.eam.fs"
lmpcmds = [
    "pair_style eam/fs",
    f"pair_coeff * * {eam_path} W",
]
calc = LAMMPSlib(lmpcmds=lmpcmds, atom_types={"W": 1}, keep_alive=True,
                 comm=worker_comm)

# ---------------------------------------------------------------------------
# Build PAFIMPI — wraps one PAFI instance per rank
# ---------------------------------------------------------------------------
pafi = PAFIMPI(images, calc, worker_comm=worker_comm, inter_comm=inter_comm)

if os.path.exists(json_dump_file):
  pafi.read_json(json_dump_file)

# Plot the NEB path once (rank 0 only)
pafi.plot_neb("w-eam-neb.png")

temperatures = [50.,500.0]
r_values = np.linspace(0.0, 0.5, 6)

for T in temperatures:
    for r in r_values:
        pafi.run(r, T=T, nsteps=30, thermsteps=30, minsteps=5, verbose=(rank == 0))
        # rank 0 now holds all N_ranks samples; write and reset on every rank
        if rank == 0:
            r_arr = np.array(pafi.data["r"])
            T_arr = np.array(pafi.data["target_T"])
            dF_arr = np.array(pafi.data["dF"])
            mask = np.isclose(r_arr,r) * np.isclose(T_arr,T) 
            print(f"  r={r:.2f}  T={T:.0f}K  n_samples={mask.sum()}  <dF>={dF_arr[mask].mean():.4f} eV")
            pafi.write_json(json_dump_file) # update as we go...
