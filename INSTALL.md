# Installation

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

## With conda
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

## From source
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

### C++ Implementation
See <a href="https://github.com/tomswinburne/pafi/tree/cpp-2023">here</a> (`cpp-2023` branch) for an older C++ implementation.
