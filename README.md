# IN BETA : FRAGILE WITH MPI/PYTHON ARCHITECTURE
         _______      _______      _______     _________
        (  ____ )    (  ___  )    (  ____ \    \__   __/
        | (    )|    | (   ) |    | (    \/       ) (
        | (____)|    | (___) |    | (__           | |
        |  _____)    |  ___  |    |  __)          | |
        | (          | (   ) |    | (             | |
        | )          | )   ( |    | )          ___) (___
        |/           |/     \|    |/           \_______/
        Projected    Average      Force        Integrator


:copyright: TD Swinburne and M-C Marinica 2018



Beta version of code used in [this paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.135503)
> *Unsupervised Calculation of Free Energy Barriers in Large Crystalline Systems*   
> T.D. Swinburne and M.-C. Marinica, Physical Review Letters 120 (13), 135503, 2018



## Patching and compiling LAMMPS

1. [Download LAMMPS](http://lammps.sandia.gov/download.html) and read the [installation instructions](http://lammps.sandia.gov/doc/Section_start.html)

2. In root directory of a clean LAMMPS repository (no packages) run
```
patch -p0 < /path/to/PAFI/lammps13Mar18_PAFI.patch
```
3. Install any packages you desire
```
make yes-package_name
```
3. Compile as shared library and copy to PAFI repo
```
   make mpi mode=shlib
   cp liblammps_mpi.so path/to/PAFI/pafilib/liblammps_mpihp.so
```



## Python Requirements
- Python >= 2.7  (< 3.0) 
- [mpi4py](mpi4py.scipy.org) >= 2.0.0
- numpy >= 1.13.0
- scipy >= 0.18.0



## Calculation of free energy barrier between states

1. First set up a LAMMPS neb calculation as described [here](http://lammps.sandia.gov/doc/neb.html)

2. In the LAMMPS script, use the "write_data" command to dump all the NEB knots, i.e.
```
variable u uloop N_NEB_IMAGES

neb etol ftol N1 N2 Nevery file-style arg keyword

write_data neb_knot_file.$u
```
3. Move the knots to folder and include a LAMMPS input script to load the first knot, as shown in the examples

4. Fill in the input file as shown in the example entry, including a wildcard command for the knots

5. Run PAFI as e.g.
```
mpirun -np NPROCS python ./lmp_mpi_hp.py -i in_file -t TEMPERATURE
```

## Output

``` 
dump_folder/free_energy_profile
```
The integrated force using the naive projection (dU/dr) and the PAFI projection (shown to be true free energy)
