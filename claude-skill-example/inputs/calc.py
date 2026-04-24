"""User-supplied calculator loader.

The `pafi-setup` skill imports `make_calculator` from this file and passes
its own `worker_comm` as the `comm=` argument. The skill does NOT rewrite
this code — if the potential / pair_style is wrong, the PAFI run is wrong.
"""

from ase.calculators.lammpslib import LAMMPSlib

# Path to the EAM potential — resolved relative to wherever run_pafi.py is
# launched from. Adjust to an absolute path if the submission script cd's
# into a different directory.
EAM_PATH = "../../neb_paths/W.eam.fs"


def make_calculator(comm=None):
    lmpcmds = [
        "pair_style eam/fs",
        f"pair_coeff * * {EAM_PATH} W",
    ]
    return LAMMPSlib(
        lmpcmds=lmpcmds,
        atom_types={"W": 1},
        keep_alive=True,
        comm=comm,
    )
