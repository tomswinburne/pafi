<img src="https://raw.githubusercontent.com/tomswinburne/pafi/refs/heads/master/doc/pafi_title.png" width=500></img>
<h2> PAFI: Evaluation of free energy barriers beyond Harmonic TST</h2>
<h4 align="center">Swinburne and Marinica, Phys. Rev. Lett 2018 (<a href="#citation">bibtex citation</a>).</h4>
PAFI performs constrained sampling on <a href="https://docs.lammps.org/fix_neb.html" target="_new">NEB</a> hyperplanes in <a href="https://docs.lammps.org" target="_new">LAMMPS</a>, 
analytically reformulating an exact expression for the free energy gradient used in the
<a href="https://pubs.acs.org/doi/10.1021/jp506633n" target="_new">Adaptive Biasing Force</a> method.
This allows calculation of free energy barriers even when the minimum energy path (MEP)
is not aligned with the minimum free energy path (MFEP). PAFI thus performs
<a href="https://en.wikipedia.org/wiki/Stratified_sampling" target="_new">stratified sampling</a> of configuration space for a particular metastable pathway, with the usual reductions in variance.
<h3 align="center">
<a href="#quick-start">Conda installation</a>
| <a href="#running-pafi">Running PAFI</a>
| <a href="#plotting-results">Plotting Results</a>
| <a href="examples/README.md">Examples</a>
| <a href="#full-installation">Full installation</a>
| <a href="#hints-and-tips">Hints and tips</a>
| <a href="#citation">Citation</a>
</h3>
</br>

## Quick Usage
If you can already run
```python
from mpi4py import MPI
from lammps import lammps
lmp = lammps()
lmp.close()
```
Then you can install PAFI:
```bash
pip install pafi
pafi-check-deps
pafi-run-tests
```

Otherwise, we can use `conda-lammps`, but **best for local testing only!** 
```bash
conda config --add channels conda-forge # add conda-forge channel
conda create -n pafi-env python=3.10 
conda activate pafi-env # activate virtual env
conda install numpy scipy pandas # install requirements (can use pip)
conda install  mpi4py lammps # conda-lammps has no MPI: one core/worker!
pip install pafi
pafi-check-deps # ensure lammps can be loaded
pafi-run-tests # run tests
```
See <a href="#full-installation">full installation</a> for optimal use on HPC. 

## Running PAFI
PAFI requires that you have already made some NEB calculation with some potential. You can then run
```shell
  mpirun -np ${NUM_PROCS} python simple_pafi_run.py
```
`simple_pafi_run.py`:
```python
  from mpi4py import MPI
  from pafi import PAFIManager, PAFIParser

  parameters = PAFIParser()
  parameters.set_potential(["neb_pathway/Fe.eam.fs"], # List of potential files
                                            "eam/fs", # LAMMPS pair style
                                             ["Fe"])  # List of Elements
  
  parameters.set_pathway("neb_pathway/image_*.dat") # NEB pathway of LAMMPS dat files
  

  parameters.axes["Temperature"] = [100.*i for i in range(7)] # temperature range
  parameters.set("CoresPerWorker",1) # Will force to 1 for conda installation of lammps
  
  # Typical production values
  parameters.set("SampleSteps",2000) # Sampling steps per worker
  parameters.set("ThermSteps",1000) # Thermalization steps per worker
  parameters.set("ThermWindow",500) # Averaging window to check temperature
  parameters.set("nRepeats",1) # Number of times to repeat process (if CPU limited)

  # Establish PAFIManager
  manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)
  manager.run()
  manager.close()
```
This will return a `*.csv` file with data readable by `pandas`. Can easily collate multiple runs. 

## Plotting Results
We can load in multiple `csv` files from some run then plot in python:
```python 
import matplotlib.pyplot as plt
from glob import glob
from pafi import ResultsProcessor

p = ResultsProcessor(   data_path=glob("dumps/pafi_data_*.csv"),
                        xml_path='dumps/config_0.xml', integrate=True)

plotting_data,x_key,y_key = p.plotting_data()


# Plotting
fig = plt.figure(figsize=(5,3))
ax = fig.add_subplot(111)
ax.set_title("Short (10ps) test for vacancy in W, EAM potential")


for i,row in enumerate(plotting_data):
    # Plotting data    
    x,y,e = row[x_key], row[y_key], row[y_key+"_err"]
    T = int(row['Temperature'])
    label = r"$\Delta\mathcal{F}$=%2.2f±%2.2f eV, T=%dK"% (y.max(),e.max(),T)
    
    # make plots
    ax.fill_between(x,y-e,y+e,facecolor=f'C{i}',alpha=0.5)
    ax.plot(x,y,f'C{i}-',lw=2,label=label)

# save
ax.legend(loc='best')
ax.set_xlabel("Reaction coordinate")
ax.set_ylabel("Free energy barrier (eV)")
```
<img src="https://raw.githubusercontent.com/tomswinburne/pafi/refs/heads/master/doc/test_output.png" width=300></img>


See the [examples](examples/README.md) and <a href="#hints-and-tips">hints and tips</a> for more information

## Full installation
See <a href="https://github.com/tomswinburne/pafi/tree/cpp-2023">here</a> (`cpp-2023` branch) for an older C++ implementation.

PAFI uses `mpi4py`, `numpy`, `scipy`, `pandas` and <b><a href="https://docs.lammps.org/Python_head.html" target="_new">LAMMPS-Python</a></b> with at least `MANYBODY` and `ML-SNAP`
If you have cmake and mpi installed:
```shell
export PREFIX=${HOME}/.local # example
export PYTHON=`which python` # to ensure same distribution
export MPICC=`which mpicc` # for mpi4py, your C++ MPI compiler (e.g. mpicc / mpiicc for intel)

# extract typical install location PLEASE CHECK THIS ON YOUR MACHINE!
# (see below for why this hack can be useful)
PYTHON_VERSION=`python --version | cut -f2 -d" " | cut -f2 -d"."`
export INSTALL_LOCATION=${PREFIX}/lib/python3.${PYTHON_VERSION}/site-packages

# get LAMMPS and PAFI source
git clone https://github.com/lammps/lammps.git
git clone https://github.com/tomswinburne/pafi.git

# install python packages
${PYTHON} -m pip install mpi4py numpy pandas

# LAMMPS build 
cd /path/to/lammps
mkdir build
cd build
cmake -C ../../pafi/doc/lammps_options.cmake ../cmake
make -j

# LAMMPS python install: 
# whilst official command is 'make install python', can have env clashes
# instead, we do it "by hand":
cd ../python # within LAMMPS repository
${PYTHON} -m pip install -U .

# manually provide binary for LAMMPS package
cp ../build/liblammps.so ${INSTALL_LOCATION}/lammps

# Install and test PAFI
cd /path/to/pafi
pip install .
python pafi/run_tests.py
```

## Hints and Tips
- See <a href="http://lammps.sandia.gov/doc/neb.html" target="_new">LAMMPS NEB</a> for making a pathway

- New `equal` style gives optimal spacing for force integration- e.g. `fix neb all neb 1.0 parallel equal`

- Modify one of `examples/configuration_files/*_REAL.xml` to load in your pathway and potential:
```python
  from mpi4py import MPI
  from pafi import PAFIManager
  manager = PAFIManager(MPI.COMM_WORLD,"/path/to/config.xml")
  manager.run()
  manager.close()
  ```

- See the [tutorial](TUTORIAL.md) for information on the `pafi-path-test` routine

- In general, we want a reference pathway with dense discretisation where energy gradients are large

- The current non-smoothed spline implementation can oscillate between very similar image configurations, as a result, there should be non-negligible displacement between images

- If your path isn't loading, try setting `LogLammps=1` in `config.xml` to check for bugs in `log.lammps`

- If `SampleSteps` is too large workers will make thermally activated "jumps" to nearby paths in the hyperplane. This will return a warning message `Reference path too unstable for sampling.`
 and increase error. If this happens, decrease `SampleSteps` and increase `nRepeats`

- When running on `NPROCS` cores, we require `NPROCS%CoresPerWorker==0`, so we have an integer number of workers

- The total number of force calls *per worker* is `nPlanes * (ThermSteps+SampleSteps) * nRepeats`, spatially parallelised by LAMMPS across `CoresPerWorker` cores for each worker.

- Each PAFI worker runs at the same speed as LAMMPS. Increasing `CoresPerWorker` will typically decrease execution time but also reduce `nWorkers` and increase error, as we have less samples.

- If you are core-limited, the `nRepeats` option forces workers to perform multiple independent sampling runs on each plane. For example, with all other parameters fixed, running on 32 cores with `nRepeats=3` is equivalent to running on 3*32=96 cores with  `nRepeats=1`, but the latter will finish in a third of the time.


## Citation
For more details please see <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.135503" target="_new">our paper</a>, citation:
```bibtex
@article{PhysRevLett.120.135503,
  title = {Unsupervised Calculation of Free Energy Barriers in Large Crystalline Systems},
  author = {Swinburne, Thomas D. and Marinica, Mihai-Cosmin},
  journal = {Phys. Rev. Lett.},
  volume = {120},
  issue = {13},
  pages = {135503},
  numpages = {6},
  year = {2018},
  month = {Mar},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.120.135503},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.120.135503}
}
```


## TODO
1. Restart files from pathway deviations
2. Smoothed spline interpolation for more general reference pathways

