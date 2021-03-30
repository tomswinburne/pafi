         _______      _______      _______     _________
        (  ____ )    (  ___  )    (  ____ \    \__   __/
        | (    )|    | (   ) |    | (    \/       ) (
        | (____)|    | (___) |    | (__           | |
        |  _____)    |  ___  |    |  __)          | |
        | (          | (   ) |    | (             | |
        | )          | )   ( |    | )          ___) (___
        |/           |/     \|    |/           \_______/
        Projected    Average      Force        Integrator

## [Main Page](README.md)

## [Installation](INSTALL.md)

## Testing the pathway

0. Steps 1,2 can be skipped if using the NEB calculation (SIA in EAM-Fe) in the examples folder

1. First set up a LAMMPS neb calculation as described [here](http://lammps.sandia.gov/doc/neb.html)

2. In the LAMMPS script, use the `write_data` command to dump all the NEB knots, i.e.
```
variable u equal part # partition for NEB image

neb etol ftol N1 N2 Nevery file-style arg keyword

write_data neb_knot_file.$u
```

3. Configure `config.xml` to load in your NEB images with the correct potential

4. Run the test routine `pafi-path-test`, which runs PAFI with a single worker,
at zero temperature, for one step per plane. Run with `CoresPerWorker` procs, e.g.
```bash
mkdir -p dumps
mpirun -np 2 ./pafi-path-test
```
where the first line ensures your dump folder (here the default value) actually exists.

5. The output of `pafi-path-test` checks the discretisation and force integration.
If warnings are raised (e.g. too few images, force integration error)
follow the suggestions to reduce error.

6. If your configuration isn't loading, try setting `LogLammps==1` in `config.xml` to find bugs


## Production Calculation

1. See tips on the [main page](README.md) for recommended values for `ThermSteps`, `SampleSteps` etc.

2. PAFI will try to write to the directory as specified in `DumpFolder` in config.xml. Each dump file has a suffix `_T_n`, where `T` is the temperature and `n` is the smallest integer that does not overwrite previous files.

3. For each temperature there will be a files `dev_r_T_n.dat` with the ensemble average and variance pathway deviation from each hyperplane and a file `free_energy_profile_T_` that has the integrated FEP.

4. See `error_analysis.pdf' for an expanation of the error bars used in PAFI and `example/sample_plot.py' for a simple plotting example

## Manual use (not recommended)
1. Compile the `pafi-lammps-preparation` binary
```bash
cd build
make pafi-lammps-path
```
2. Configure the configuration xml file, to specify the path then run
```bash
mkdir -p dumps
mpirun -np 1 ./pafi-lammps-preparation
```

3. This will make a set of files `dumps/pafipath.*.data` which can be run using the `PAFI` example script in the `examples/USER/misc/pafi` directory of the `LAMMPS` repository
