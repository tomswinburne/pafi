from mpi4py import MPI
from pafi import PAFIManager,PAFIParser

"""
    Overwrite parameters within python script
"""

# load in default parameters
parameters = PAFIParser()
# set path wildcard, assuming some_file_name_INTEGER.dat format
parameters.set_pathway("systems/EAM-SIA-Fe/image_*.dat")
# set interatomic potential
parameters.set_potential("systems/EAM-SIA-Fe/Fe.eam.fs")
# restrict to zero temperature
parameters.axes["Temperature"] = [100.]
parameters.set("nRepeats",2)
parameters.set("SampleSteps",10)
parameters.set("ThermSteps",10)
parameters.set("ThermWindow",10)

manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)
manager.run()
manager.close()