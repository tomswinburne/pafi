from mpi4py import MPI
from pafi import PAFIManager,PAFIParser

# Overwrite parameters within default python script (assumes eam/fs)
parameters = PAFIParser()
parameters.set_pathway("systems/EAM-VAC-W/image_*.dat")
parameters.set_potential("systems/EAM-VAC-W/W.eam.fs",type="eam/fs")
parameters.set_species("W")

# TEST VALUES
parameters.axes["Temperature"] = [100.]
parameters.set("nRepeats",1)
parameters.set("SampleSteps",10)
parameters.set("ThermSteps",10)
parameters.set("ThermWindow",10)

manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)
manager.run()
manager.close()
