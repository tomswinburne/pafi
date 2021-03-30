         _______      _______      _______     _________
        (  ____ )    (  ___  )    (  ____ \    \__   __/
        | (    )|    | (   ) |    | (    \/       ) (
        | (____)|    | (___) |    | (__           | |
        |  _____)    |  ___  |    |  __)          | |
        | (          | (   ) |    | (             | |
        | )          | )   ( |    | )          ___) (___
        |/           |/     \|    |/           \_______/
        Projected    Average      Force        Integrator

# Evaluate free energy barriers beyond the harmonic approximation

## "Dropping the H in HTST"

v0.9 :copyright: TD Swinburne and M-C Marinica 2020 MIT License

swinburne at cinam.univ-mrs.fr

Beta version of code used in [this paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.135503)
> *Unsupervised Calculation of Free Energy Barriers in Large Crystalline Systems*   
> T.D. Swinburne and M.-C. Marinica, Physical Review Letters 120 (13), 135503, 2018

Please cite the above when publishing results using PAFI

This repository includes the [RapidXML](http://http://rapidxml.sourceforge.net) library for input parsing

## [Installation Instructions](INSTALL.md)

## [Getting Started](TUTORIAL.md)

## Notes on reference pathway

- See the [tutorial](TUTORIAL.md) for information on the `pafi-path-test` routine

- In general, we want a reference pathway with dense discretisation where energy gradients are large, and vice versa

- The current non-smoothed spline implementation (see TODO) can oscillate between two very similar image configurations, as a result, there should be non-negligible displacement between images

- If your path isn't loading, set `LogLammps=1` in `config.xml` to check for bugs in `log.lammps`

## Notes on choosing parameters

- If `SampleSteps` is too large workers will make thermally activated "jumps" to nearby paths in the hyperplane. This will return a warning message `Reference path too unstable for sampling.`
 and increase error. If this happens, decrease `SampleSteps` and increase `nRepeats`

- When running on `NPROCS` cores, we require `NPROCS%CoresPerWorker==0`, so we have an integer number of workers

- The total number of force calls *per worker* is `nPlanes * (ThermSteps+SampleSteps) * nRepeats`, spatially parallelised by LAMMPS across `CoresPerWorker` cores for each worker.

- Each PAFI worker runs at the same speed as LAMMPS. Increasing `CoresPerWorker` will typically decrease execution time but also reduce `nWorkers` and increase error, as we have less samples.

- If you are core-limited, the `nRepeats` option forces workers to perform multiple independent sampling runs on each plane. With all other parameters fixed, `nRepeats=3, NPROCS=32` is equivalent to `nRepeats=1, NPROCS=96`

## TODO
1. Restart files from pathway deviations
2. Smoothed spline interpolation for more general reference pathways
