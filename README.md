<img src="doc/pafi_title.png" width=500></img>
<h2> PAFI: Evaluation of free energy barriers beyond Harmonic TST</h2>
<h3 align='center'>
Python implementation (2024), (c) Swinburne & Marinica 2024<br><br>
<a href="https://github.com/tomswinburne/pafi/tree/cpp-2023">C++ implementation (2023) available here</a>
<h3 align="center">If using PAFI, please cite (<a href="#citation">bibtex</a>). Details in 
<a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.135503" target="_new">Swinburne & Marinica PRL 2018</a>.
</h3>
PAFI performs constrained sampling on <a href="https://docs.lammps.org/fix_neb.html" target="_new">NEB</a> hyperplanes in 
<a href="https://docs.lammps.org" target="_new">LAMMPS</a>, 
analytically reformulating an exact expression for the free energy gradient used in the
<a href="https://pubs.acs.org/doi/10.1021/jp506633n" target="_new">Adaptive Biasing Force</a> method.
This allows calculation of free energy barriers even when the minimum energy path (MEP)
is not aligned with the minimum free energy path (MFEP). PAFI thus performs
<a href="https://en.wikipedia.org/wiki/Stratified_sampling" target="_new">stratified sampling</a> of configuration 
space for a particular metastable pathway, with the usual reductions in variance.
</br>
</br>


# Quick Installation
- Python 3 with `numpy`, `mpi4py` and `tqdm` (optional)
- <b><a href="https://docs.lammps.org/Python_head.html" target="_new">LAMMPS-Python</a></b> with at least `MANYBODY` and `ML-SNAP`
- You should be able to run the following python code:
	```python
	from mpi4py import MPI
	from lammps import lammps
	comm = MPI.COMM_WORLD
	lmp = lammps(comm=comm)
	lmp.close()
	```
- Local install with `pip` and testing:
	```bash
	cd /path/to/DeFAD # this repo
	# conda activate defad_env # or similar, if using venv
	python -m pip install -e . # local install
  cd testing # test suite
  python run_tests.py -v
	```

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

## Full installation
If you have cmake and mpi installed:
```bash
export PREFIX=${HOME}/.local # example
export PYTHON=`which python` # to ensure same distribution
export MPICC=`which mpicc` # for mpi4py, your C++ MPI compiler (e.g. mpicc / mpiicc for intel)

# extract typical install location PLEASE CHECK THIS ON YOUR MACHINE!
# (see below for why this hack can be useful)
PYTHON_VERSION=`python --version | cut -f2 -d" " | cut -f2 -d"."`
export INSTALL_LOCATION=${PREFIX}/lib/python3.${PYTHON_VERSION}/site-packages

# get souces
git clone https://github.com/lammps/lammps.git
git clone https://github.com/tomswinburne/pafi.git

# install python packages
${PYTHON} -m pip install mpi4py numpy tqdm

# LAMMPS build 
cd /path/to/lammps
mkdir build
cd build
cmake -C ../../pafi/cmake/lammps_options.cmake ../cmake
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
pip install -e .
cd testing # within PAFI repository
python run_tests.py

```

## [Detailed Installation Instructions](doc/INSTALL.md)
## [Getting Started Tutorial](doc/TUTORIAL.md)
## [Hints and Tips](doc/TIPS.md)

## TODO
1. Restart files from pathway deviations
2. Smoothed spline interpolation for more general reference pathways
3. More unit tests...
4. incorporate Arnauds path preparation scripts

