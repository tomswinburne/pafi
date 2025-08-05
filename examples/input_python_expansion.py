import os
from glob import glob
from mpi4py import MPI
from pafi import PAFIManager,PAFIParser

# Overwrite parameters within default python script (assumes eam/fs)
rank = MPI.COMM_WORLD.Get_rank()
parameters = PAFIParser(rank=rank)

# force field
parameters.set_potential("systems/EAM-VAC-W/W.eam.fs",pot_type="eam/fs")
parameters.set_species("W")

# set file list 
file_dir = "systems/EAM-EXPAND-W"
file_wildcard = "strain_*.dat"
# function to extract index
file_index = lambda f:int(f.split(".")[-2].split("_")[-1])
# list all files
file_list = glob(os.path.join(file_dir,file_wildcard))
# remove location from path
file_list = [os.path.basename(f) for f in file_list]
# sort according to index
file_list = sorted(file_list, key=file_index)
# add to parameters
parameters.set_pathway(file_list, directory=file_dir)

# TEST VALUES
parameters.axes["Temperature"] = [0.0,300.,1000.,2000.]
parameters.set("nRepeats",1)
parameters.set("OverDamped",0)
parameters.set("SampleSteps",2000)
parameters.set("ThermSteps",1000)
parameters.set("ThermWindow",100)
parameters.set("LogLammps",0)
parameters.set("PostDump",0)

manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)
manager.run()
manager.close()
