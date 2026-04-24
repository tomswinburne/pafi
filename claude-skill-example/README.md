# `pafi-setup` skill — worked example

This folder shows exactly what the `pafi-setup` Claude Code skill consumes
and what it produces. Nothing here is executed by the skill — it's a
reference so you can see the shape of the inputs and the shape of the
generated outputs before using it on your own system.

## Layout

```
claude-skill-example/
├── inputs/          # what the user places in their working folder
│   ├── calc.py          # calculator loader  (LAMMPSlib / W EAM)
│   ├── neb.traj         # NEB pathway        (symlink to repo's w-vac-path.traj)
│   └── submit.slurm     # scheduler template (SLURM, partially filled)
└── outputs/         # what the skill writes into the user's folder
    ├── run_pafi.py      # generated PAFI driver
    └── submit_pafi.sh   # generated, filled-in submission script
```

## The scenario

A user has a tungsten vacancy migration pathway (the `w-vac-path.traj` that
ships with this repo) and wants to run PAFI on a SLURM cluster using a
LAMMPS EAM potential. They drop three files into a working folder and ask
Claude Code:

> set up a PAFI run in this folder — calc.py has the calculator, neb.traj is
> the NEB, use submit.slurm as the template, grid r=0..1 in 11 points at
> T=300K, smoke-test sizes

## Inputs (`inputs/`)

- **`calc.py`** — exposes `make_calculator(comm=None)` returning an ASE
  `LAMMPSlib` instance. The skill imports this symbol; it does not rewrite
  calculator construction.
- **`neb.traj`** — a symlink to `../../neb_paths/w-vac-path.traj`. Any
  `ase.io.read(..., index=":")`-compatible file works (traj, xyz, extxyz).
- **`submit.slurm`** — a realistic SLURM header with `#SBATCH` directives,
  a module load block, and a `<LAUNCH>` placeholder on the launch line. The
  skill preserves the user's directives and only fills the placeholder.

## Outputs (`outputs/`)

- **`run_pafi.py`** — MPI driver. Imports `make_calculator` from `calc.py`,
  loads `neb.traj`, wires `PAFIMPI` with matching `worker_comm`, and runs
  the requested `(r, T)` grid with JSON restart.
- **`submit_pafi.sh`** — the user's `submit.slurm` with its directives kept
  verbatim, placeholders resolved, and the launch line replaced by
  `mpirun -n $SLURM_NTASKS python run_pafi.py`.

The skill would also print a launch command at the end:

```
sbatch submit_pafi.sh
```

It never submits the job itself.
