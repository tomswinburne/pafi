         _______      _______      _______     _________
        (  ____ )    (  ___  )    (  ____ \    \__   __/
        | (    )|    | (   ) |    | (    \/       ) (
        | (____)|    | (___) |    | (__           | |
        |  _____)    |  ___  |    |  __)          | |
        | (          | (   ) |    | (             | |
        | )          | )   ( |    | )          ___) (___
        |/           |/     \|    |/           \_______/
        Projected    Average      Force        Integrator

# MD evaluation of free energy barriers beyond HTST
v0.9 :copyright: TD Swinburne and M-C Marinica 2020 MIT License

thomas dot swinburne at cnrs dot fr

:rotating_light: New `equal` NEB style for `LAMMPS` on [github](https://github.com/lammps/lammps.git), ideal to produce pathways for PAFI [documentation](https://github.com/lammps/lammps/blob/develop/doc/src/fix_neb.rst) :rotating_light:

Using PAFI? Please cite [this paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.135503)
> *Unsupervised Calculation of Free Energy Barriers in Large Crystalline Systems*   
> T.D. Swinburne and M.-C. Marinica
> Physical Review Letters 120 (13), 135503, 2018

Applications:
> Namakian *et al.* Comp. Mat. Sci. 2023 [link](https://doi.org/10.1016/j.commatsci.2022.111971)

> Sato *et al.* Mat. Res. Lett.  2021 [link](https://doi.org/10.1080/21663831.2021.1875079)


## [Installation Instructions](INSTALL.md)

## [Getting Started Tutorial](TUTORIAL.md)

## General Tips

- See the [tutorial](TUTORIAL.md) for information on the `pafi-path-test` routine

- In general, we want a reference pathway with dense discretisation where energy gradients are large

- The current non-smoothed spline implementation can oscillate between very similar image configurations, as a result, there should be non-negligible displacement between images

- If your path isn't loading, try setting `LogLammps=1` in `config.xml` to check for bugs in `log.lammps`

## Notes on choosing parameters

- The `RealMEPDist` parameter determines how the "reaction coordinate" (r.c.) is calculated. For  `RealMEPDist=0` the r.c. is simply `HyperplaneIndex/(nPlanes-1)`. `RealMEPDist=1` the r.c. is the total real space distance from the first configuration, normalized to one, whilst for `RealMEPDist=2`, the r.c. is the average distance from the first and to the last configuration, normalized to one. For best comparision to a NEB calculation, we recommend `RealMEPDist=1`.

- If `SampleSteps` is too large workers will make thermally activated "jumps" to nearby paths in the hyperplane. This will return a warning message `Reference path too unstable for sampling.`
 and increase error. If this happens, decrease `SampleSteps` and increase `nRepeats`

- When running on `NPROCS` cores, we require `NPROCS%CoresPerWorker==0`, so we have an integer number of workers

- The total number of force calls *per worker* is `nPlanes * (ThermSteps+SampleSteps) * nRepeats`, spatially parallelised by LAMMPS across `CoresPerWorker` cores for each worker.

- Each PAFI worker runs at the same speed as LAMMPS. Increasing `CoresPerWorker` will typically decrease execution time but also reduce `nWorkers` and increase error, as we have less samples.

- If you are core-limited, the `nRepeats` option forces workers to perform multiple independent sampling runs on each plane. For example, with all other parameters fixed, running on 32 cores with `nRepeats=3` is equivalent to running on 3*32=96 cores with  `nRepeats=1`, but the latter will finish in a third of the time.



## External Libraries
- [LAMMPS](https://lammps.sandia.gov) MD code
- [RapidXML](https://rapidxml.sourceforge.net) for reading `xml` files
- This nice [library](https://github.com/ttk592/spline) for spline interpolation

## TODO
1. Restart files from pathway deviations
2. Smoothed spline interpolation for more general reference pathways
