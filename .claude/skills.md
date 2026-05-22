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
pip install --force-reinstall --no-deps .   # or: pip install -e .
```

from the repo root. Note this is **not** an editable install unless you pass `-e`; with `--force-reinstall --no-deps .` the worker Python processes spawned by `mpirun` will resolve `pafi` from site-packages, not from the working tree, so any subsequent edits won't take effect without reinstalling. If you're iterating on source, use `pip install -e .` (or set `PYTHONPATH` to the repo root).

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

## Generating new driver scripts

When a user says something like *"make me a PAFI example for my <potential> on <system>"* (or any variant — "write a driver", "generate a script", "set this up for ACE NiAl", …), follow the routine here.

### Info to gather

If the user hasn't supplied it (and you can't infer from a filename or directory), ask for:

1. **Potential file(s)** — path(s) to the parameter file. The extension is usually enough to guess the pair_style; confirm.
2. **Species list, in order** — must match the `read_data` atom types (type 1 = first species, etc.) and the `pair_coeff` line.
3. **NEB pathway** — glob like `image_*.dat` plus the directory holding them. Don't guess; ask if there are non-default boundary conditions or empty `Masses` blocks.
4. **Temperature(s)** in K (one value or a list).
5. **Quality preset** — smoke / quick / production (see below).
6. **CoresPerWorker** — defaults to 1; set higher for big systems or HPC.

### Pair-style routing

Decide which template to use based on pair_style. **Simple path**: PAFI's default `Input` script does the right thing — just call `set_potential()` + `set_species()` and stop. **Custom path**: write your own `Input` and `PreRun` scripts; `set_potential()` becomes optional (used only for the `%Potential%` substitution token).

| pair_style | File ext. | `pair_coeff` shape | Template |
|---|---|---|---|
| `eam/fs`, `eam/alloy` | `.eam.fs`, `.eam.alloy` | `* * <file> <species…>` | **Simple** |
| `eam` (old Funcfl) | `.eam` | `<i> <j> <file>` per pair | Custom |
| `snap` | `.snapcoeff` + `.snapparam` | `* * <coeff> <param> <species…>` | Custom |
| `mlip` / `pace` / `grace/fs` | `.yaml`, `.ace` | `* * <file> <species…>` | Custom (non-default pair_style name) |
| `hybrid` / `hybrid/scaled` | mixed | `* * <sub_style> <args>` per pair | Custom |
| `lj/cut`, `morse`, etc. | none | `<i> <j> <eps> <sig> …` | Custom (no file) |

Rule of thumb: if `pair_coeff` is exactly `* * <one_file> <species>` *and* the pair_style is one of `eam/fs`/`eam/alloy`, use the Simple template. Anything else → Custom.

### Template A — Simple (eam/fs, eam/alloy)

```python
from mpi4py import MPI
from pafi import PAFIManager, PAFIParser

rank = MPI.COMM_WORLD.Get_rank()
parameters = PAFIParser(rank=rank)

parameters.set_pathway("image_*.dat", directory="./neb_path")
parameters.set_potential("./pots/Fe.eam.fs", pot_type="eam/fs")
parameters.set_species(["Fe"])              # order = read_data atom types

parameters.axes["Temperature"] = [300.0]
parameters.set("CoresPerWorker", 1)
parameters.set("nRepeats", 1)
parameters.set("OverDamped", 0)
parameters.set("SampleSteps", 2000)
parameters.set("ThermSteps", 2000)
parameters.set("ThermWindow", 100)

manager = PAFIManager(MPI.COMM_WORLD, parameters=parameters)
manager.run()
manager.close()
```

### Template B — Custom Input + PreRun (ACE / GRACE / SNAP / hybrid / unusual pair_style)

Use this when the default `Input` script's `pair_style %PotentialType%` line wouldn't be valid LAMMPS (e.g. `grace/fs` needs explicit args, `snap` needs two files, `hybrid` needs sub-style specs).

```python
from mpi4py import MPI
from pafi import PAFIManager, PAFIParser

rank = MPI.COMM_WORLD.Get_rank()
parameters = PAFIParser(rank=rank)

parameters.set_pathway("image_*.dat", directory="./dat_files")
parameters.set_potential("FS_model.yaml", pot_type="grace/fs")  # optional; enables %Potential%
parameters.set_species(["Ni", "Al"])         # order matters: matches mass/pair_coeff below

parameters.axes["Temperature"] = [300.0]
parameters.set("CoresPerWorker", 64)
parameters.set("nRepeats", 1)
parameters.set("OverDamped", 1)
parameters.set("SampleSteps", 2000)
parameters.set("ThermSteps", 2000)
parameters.set("ThermWindow", 100)
parameters.set("LogLammps", 1)

parameters.set_script("Input", """
    units           metal
    dimension       3
    boundary        p p s
    atom_style      atomic
    atom_modify     map array sort 0 0.0
    newton          on
    neigh_modify    every 2 delay 10 check yes page 1000000 one 100000
    read_data       %FirstPathConfiguration%
    mass 1 58.69    # Ni
    mass 2 26.98    # Al
""")

parameters.set_script("PreRun", """
    pair_style      grace/fs
    pair_coeff      * * FS_model.yaml Ni Al
""")

manager = PAFIManager(MPI.COMM_WORLD, parameters=parameters)
manager.run()
manager.close()
```

Notes for Template B:
- **Why `pair_style`/`pair_coeff` go in `PreRun`, not `Input`**: PAFI re-runs `PreRun` at each hyperplane; if you ever need per-hyperplane pair changes (thermal expansion, switching potentials, etc.) you do it here. For a static potential it works equally in either, but matching the existing custom example keeps things consistent.
- **`mass` lines**: needed when the LAMMPS `.dat` file's `Masses` block is missing or zero. Check `head -25 image_0.dat`.
- **`boundary`**: default `p p p`; use `p p s` (slab) or `f f f` (cluster) when the system demands it. Ask the user if you can't tell.
- **Species order**: `set_species(["Ni","Al"])` must match `mass 1 Ni / mass 2 Al` and `pair_coeff * * <file> Ni Al`. Off-by-one here causes silent wrong physics.

### Available substitution tokens

Within `set_script()` bodies, the parser replaces `%Token%` with current values. Available everywhere:

| Token | Source | When |
|---|---|---|
| `%FirstPathConfiguration%` | `PathwayConfigurations[0]` | Always |
| `%Potential%` | `set_potential(path)` | If `set_potential` was called |
| `%PotentialType%` | `set_potential(..., pot_type=...)` | If `set_potential` was called |
| `%Species%` | `set_species([...])` joined by spaces | If `set_species` was called |

Additional tokens in per-hyperplane scripts (`PreRun`, `PreTherm`, `PostTherm`, `PostRun`):

| Token | Meaning |
|---|---|
| `%Temperature%` | Current temperature (K) for this hyperplane |
| `%ReactionCoordinate%` | Current reaction coordinate ∈ [0, 1] |
| `%SampleSteps%`, `%ThermSteps%`, `%ThermWindow%` | Override to 1 automatically at T=0 |
| `%Repeat%` | 1-based repeat index (if `nRepeats > 1`) |
| Anything else in `parameters.axes` | e.g. add `parameters.axes["Strain"] = [...]` and `%Strain%` becomes available |

### Quality presets

Pick one and fill it in. Don't ship "production" defaults silently — ask first.

```python
# smoke  — single-T, runs in <1 min, just checks the script doesn't crash
parameters.axes["Temperature"] = [0.0]
parameters.set("SampleSteps", 10); parameters.set("ThermSteps", 10)
parameters.set("ThermWindow", 10); parameters.set("nRepeats", 1)

# quick  — gives a noisy barrier estimate, useful for path debugging
parameters.set("SampleSteps", 500);  parameters.set("ThermSteps", 500)
parameters.set("ThermWindow", 100); parameters.set("nRepeats", 1)

# production — what *_REAL.xml uses, minutes-to-hours per T
parameters.set("SampleSteps", 5000); parameters.set("ThermSteps", 2000)
parameters.set("ThermWindow", 500); parameters.set("nRepeats", 3)
```

### Checklist before handing the script back

- `set_species(...)` order matches the `mass N <element>` lines and the trailing species on `pair_coeff`.
- `%FirstPathConfiguration%` appears in `Input` (the parser substitutes the first image; PAFI replays `read_data` for each hyperplane internally).
- If pair_style isn't `eam/fs`/`eam/alloy`, you supplied a custom `Input` and `PreRun`.
- For non-default boundary conditions, you set `boundary` explicitly in `Input` (LAMMPS default is `p p p`).
- Run command stated: `mpirun -np $(N) python <script>.py` from the directory containing the data files (or with absolute paths in `set_pathway`).
- `CoresPerWorker` divides total MPI ranks evenly (`nWorkers = NPROCS // CoresPerWorker`).

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
| `CoresPerWorker` | 1 | LAMMPS ranks per PAFI worker. `nWorkers = NPROCS / CoresPerWorker` and the division must be exact. |
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
