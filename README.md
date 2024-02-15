# PAFI Examples

More detailed / involved examples will be added soon !

## Systems
We have two point defects relaxed with an EAM potential
- `systems/EAM-SIA-Fe` : <br>
    dumbell SIA with Marinica07 EAM Fe
- `systems/EAM-VAC-W` : <br>
    vacancy with Marinica04 EAM W

## PAFI tests 
Full unittests will be implemented in the future

! Example input files `configuration_files/PartialConfiguration_TEST.xml` and `configuration_files/CompleteConfiguration_TEST.xml` have **very short samples for testing at zero or very low temperature - DO NOT USE THESE FOR SCIENCE!** 

! Please see `configuration_files/PartialConfiguration_REAL.xml` or `configuration_files/CompleteConfiguration_REAL.xml` for realistic values

Run test using "complete" input file with four workers:
```bash
cd examples/
mpirun -np 4 python UsageExamples.py -t complete
```

Run test using python input with four workers:
```bash
cd examples/
mpirun -np 4 python UsageExamples.py -t python
```

Test postprocessing:
```bash
cd examples/
python UsageExamples.py -t integrate
```