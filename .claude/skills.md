# PAFI — chat-based documentation

Reference for running PAFI (Projected Average Force Integrator) on top of LAMMPS-Python. Captures install gotchas and the example workflow.

## Installation

### Quick pip route

```bash
pip install lammps pafi
```

This installs `lammps`, `pafi`, and pulls in `mpi4py`, `numpy`, `scipy`, `pandas`, `matplotlib`. `pip install lammps` ships a prebuilt LAMMPS wheel — no compilation needed.

**MPI runtime is still required.** `mpi4py` and the LAMMPS wheel both `dlopen` `libmpi`, and they must agree on the ABI:

- The PyPI LAMMPS wheel is linked against **MPICH** (`libmpi.so.12`).
- OpenMPI provides `libmpi.so.40` and will fail with:
  ```
  OSError: libmpi.so.12: cannot open shared object file
  ```
- On Debian/Ubuntu:
  ```bash
  apt-get install -y libmpich-dev mpich
  update-alternatives --set mpirun /usr/bin/mpirun.mpich   # if OpenMPI was already there
  ```
- If both MPICH and OpenMPI are installed, remove OpenMPI or you'll get symbol-resolution crashes when `mpi4py` and LAMMPS load different `libmpi`s.

### Verification

```bash
python -c "from mpi4py import MPI; from lammps import lammps; lmp = lammps(); lmp.close()"
pafi-check-deps
```

Both should report success. Note: `lammps()` (and therefore `pafi-check-deps`) writes a `log.lammps` to the current directory unless `-log none` is passed. PAFI worker instances are controlled separately via the `LogLammps` parameter (see below).

### From local repo

To install from a checkout instead of PyPI:

```bash
pip install --force-reinstall --no-deps .
```

from the repo root.

## Examples

All under `examples/`. Run with MPICH's `mpirun`:

```bash
cd examples
mpirun -np 2 python input_python.py
```

| Script | What it does |
|---|---|
| `input_python.py` | W vacancy, EAM, sets params in Python |
| `input_python_custom.py` | Same system with custom `Input`/`PreRun` scripts and `hybrid/scaled` pair style |
| `input_xml.py` | Reads `configuration_files/CompleteConfiguration_TEST.xml` |
| `post_processing.py` | Reads `dumps/pafi_data_*.csv`, integrates `<dF/dx>`, prints barriers |

Outputs land in `examples/dumps/`:
- `pafi_data_<N>.csv` — one per worker, per temperature
- `config_<N>.xml` — full effective configuration as written by PAFI

### Minimal Python usage

```python
from mpi4py import MPI
from pafi import PAFIManager, PAFIParser

rank = MPI.COMM_WORLD.Get_rank()
parameters = PAFIParser(rank=rank)

parameters.set_pathway("image_*.dat", directory="systems/EAM-VAC-W")
parameters.set_potential("systems/EAM-VAC-W/W.eam.fs", pot_type="eam/fs")
parameters.set_species("W")

parameters.axes["Temperature"] = [0., 1000., 2000.]
parameters.set("nRepeats", 1)
parameters.set("OverDamped", 0)
parameters.set("SampleSteps", 2000)
parameters.set("ThermSteps", 2000)
parameters.set("ThermWindow", 100)

manager = PAFIManager(MPI.COMM_WORLD, parameters=parameters)
manager.run()
manager.close()
```

### Custom LAMMPS scripts

`PAFIParser.set_script(stage, body)` injects raw LAMMPS into the worker. Stages, in order:

1. `Input` — runs once at worker start. Must `read_data %FirstPathConfiguration%`, set `pair_style`/`pair_coeff`, etc.
2. `PreRun` — before each hyperplane, before `fix pafi`
3. `PreTherm` — after `fix pafi`, before thermalisation
4. `PostTherm` — after thermalisation, before sampling
5. `PostRun` — after `unfix pafi`

The parser substitutes `%Tokens%` (e.g. `%FirstPathConfiguration%`, `%Temperature%`) into scripts.

### Post-processing

```python
from pafi import ResultsProcessor
import glob

p = ResultsProcessor(data_path=glob.glob("dumps/pafi_data_*.csv"))
_, integrated = p.integrate(
    argument="ReactionCoordinate",
    target="FreeEnergyGradient",
    remesh=10,
    return_remeshed_array=True,
)
for d in integrated:
    i = d["FreeEnergyGradient_integrated"].argmax()
    print(d["Temperature"], d["FreeEnergyGradient_integrated"][i])
```

## Key parameters

Set via `parameters.set("Name", value)` in Python or `<Name>value</Name>` inside `<Parameters>` in XML.

| Parameter | Default | Meaning |
|---|---|---|
| `Temperature` (axis) | — | List of temperatures in K. Set via `parameters.axes["Temperature"] = [...]`. |
| `nRepeats` | — | Independent samples per (T, hyperplane). |
| `SampleSteps` | — | MD steps used for the PAFI average. |
| `ThermSteps` | — | MD steps for thermalisation before sampling. |
| `ThermWindow` | — | Window for thermalisation convergence check. |
| `OverDamped` | 0 | `1` = Brownian dynamics, `0` = Langevin. |
| `LogLammps` | 0 | `1` → each worker writes `log.lammps.<worker_instance>`; `0` → `-log none`. Controls only the PAFI workers, not the standalone `lammps()` in `pafi-check-deps`. |
| `GlobalSeed` | 137 | RNG seed root. |
| `SplinePath` / `RealMEPDist` / `ReDiscretize` | 1 | Pathway interpolation toggles. |
| `LinearThermalExpansion` / `QuadraticThermalExpansion` | `[0,0,0]` | Lattice expansion vs T. |

## Test systems

Bundled in `examples/systems/`:

- `EAM-VAC-W` — vacancy in W (Marinica04 EAM), 9 images `image_0.dat`…`image_8.dat`, `W.eam.fs`.
- `EAM-SIA-Fe` — dumbbell SIA in Fe (Marinica07 EAM), 9 images, `Fe.eam.fs`.

The XML "TEST" configs use sampling values too short for science — they exist only for smoke tests. Use `*_REAL.xml` for real runs.

## Layout

```
pafi/
  parsers/BaseParser.py         # set_pathway, set_potential, set_species, set_script, read_pathway
  parsers/PAFIParser.py         # main user-facing parser
  managers/PAFIManager.py       # orchestrates workers + gatherer
  managers/BaseManager.py
  workers/LAMMPSWorker.py       # start_lammps(), run_script(), LogLammps handling at line 91
  workers/PAFIWorker.py         # PAFI-specific worker logic
examples/
  input_python.py
  input_python_custom.py
  input_xml.py
  post_processing.py
  configuration_files/*.xml
  systems/{EAM-VAC-W,EAM-SIA-Fe}/
```
