from mpi4py import MPI
from pafi import PAFIManager,PAFIParser

# Overwrite parameters within default python script (assumes eam/fs)
rank = MPI.COMM_WORLD.Get_rank()
parameters = PAFIParser(rank=rank)
parameters.set_pathway("image_*.dat",directory="systems/EAM-VAC-W")
parameters.set_potential("systems/EAM-VAC-W/W.eam.fs",pot_type="eam/fs")
parameters.set_species("W")

# TEST VALUES
parameters.axes["Temperature"] = [0.,1000.,2000.]
parameters.set("nRepeats",1)
parameters.set("OverDamped",0)
parameters.set("SampleSteps",2000)
parameters.set("ThermSteps",2000)
parameters.set("ThermWindow",100)

manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)
manager.run()
manager.close()
