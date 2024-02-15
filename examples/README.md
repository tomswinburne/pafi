# PAFI Examples

## Test Systems
- `systems/EAM-SIA-Fe` : Dumbell SIA with Marinica07 EAM Fe
- `systems/EAM-VAC-W` : Vacancy with Marinica04 EAM W

## Test scripts
- XML files can be used to overwrite default parameters, or they can be set in python

- `configuration_files/*_TEST.xml` are for quick testing- **DO NOT USE THESE FOR SCIENCE!** 

- `configurations_files/*_REAL.xml` have realistic sampling values.

- All run with e.g. `export NPROCS=4;mpirun -np ${NPROCS} python input_python.py`

- Example XML usage in `input_xml.py`:
```python
    from mpi4py import MPI
    from pafi import PAFIManager
    config = "./configuration_files/PartialConfiguration_TEST.xml"
    manager = PAFIManager(MPI.COMM_WORLD,config)
    manager.run()
    manager.close()
```

- Example python usage in `input_python.py`:
```python
    from mpi4py import MPI
    from pafi import PAFIManager
    # Overwrite parameters within python script
    parameters = PAFIParser()
    parameters.set_pathway("systems/EAM-SIA-Fe/image_*.dat")
    parameters.set_potential("systems/EAM-SIA-Fe/Fe.eam.fs")

    # TEST VALUES
    parameters.axes["Temperature"] = [100.]
    parameters.set("nRepeats",2)
    parameters.set("SampleSteps",10)
    parameters.set("ThermSteps",10)
    parameters.set("ThermWindow",10)

    manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)
    manager.run()
    manager.close()
```
## Post processing
- Can print integration results with `python post_processing.py`

- See also jupyter notebook `PlottingExamples.ipynb`
