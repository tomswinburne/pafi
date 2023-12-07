<img src="doc/pafi_title.png" width=500></img>
<h2> PAFI: MD evaluation of free energy barriers beyond HTST</h2>
v0.9 :copyright: TD Swinburne and M-C Marinica 2023 MIT License, thomas dot swinburne at cnrs.fr<br><br>
<h3>Currently in beta, stable release available at https://github.com/tomswinburne/pafi</h3>

PAFI performs constrained sampling on [NEB](https://docs.lammps.org/fix_neb.html) hyperplanes in [LAMMPS](https://docs.lammps.org), 
analytically reformulating an exact expression for the free energy gradient used in the
[Adaptive Biasing Force](https://pubs.acs.org/doi/10.1021/jp506633n) method.
This allows calculation of free energy barriers even when the minimum energy path (MEP)
is not aligned with the minimum free energy path (MFEP). PAFI thus performs
[stratified sampling](https://en.wikipedia.org/wiki/Stratified_sampling) of configuration 
space for a particular metastable pathway, with the usual reductions in variance.
For more details please see (and cite) [our paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.135503):
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

## Quick build
If you have cmake and mpi installed:
```bash
export PREFIX=${HOME}/.local # example
git clone https://github.com/lammps/lammps.git
git clone https://github.com/tomswinburne/pafi.git

# LAMMPS build
mkdir lammps/build
cd lammps/build
cmake -C ../../pafi/cmake/lammps_options.cmake ../cmake
make -j
make install # to PREFIX
make install python

# back to parent
cd ..

# PAFI build
cd pafi
pip install -e .
```

## [Detailed Installation Instructions](doc/INSTALL.md)
## [Getting Started Tutorial](doc/TUTORIAL.md)
## [Hints and Tips](doc/TIPS.md)

## External Libraries
- [LAMMPS](https://lammps.sandia.gov) MD code
- [RapidXML](https://rapidxml.sourceforge.net) for reading `xml` files
- This nice [library](https://github.com/ttk592/spline) for spline interpolation

## TODO
1. Restart files from pathway deviations
2. Smoothed spline interpolation for more general reference pathways
3. print parameter object as csv also, or json
4. have proper testing routine
5. incorporate Arnauds path preparation scripts

