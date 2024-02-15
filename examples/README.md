# PAFI Examples
Please run `testing/run_tests.py` for unit tests

## Configuration files
- `configuration_files/*_TEST.xml`: unrealistic parameters for rapid testing **DO NOT USE THESE FOR SCIENCE!** 
- `configuration_files/*_REAL.xml`: realistic parameters for most applications

## Systems
We have two point defects relaxed with an EAM potential
- `systems/EAM-SIA-Fe` : dumbell SIA with Marinica07 EAM Fe
- `systems/EAM-VAC-W` : vacancy with Marinica04 EAM W

## Example scripts
- Simple test with complete or partial configuration file:
```bash
mpirun -np 4 python standard_input.py
```

- Overwriting default values within python script:
```bash
mpirun -np 4 python python_input.py
```

- Post processing results (after running one of the above)
```bash
python post_processing.py
```

- See also jupyter notebook `PlottingExamples.ipynb`