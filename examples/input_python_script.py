from mpi4py import MPI
from pafi import PAFIManager,PAFIParser


"""
    More flexible implementation for custom pair styles

    Here, redundant example mixing the same potential twice 

"""
rank = MPI.COMM_WORLD.Get_rank()
parameters = PAFIParser(rank=rank)

parameters.set_pathway("image_*.dat",directory="systems/EAM-VAC-W")

# Currently also need to give species to PAFI here (will be made redundant) 
parameters.set_species(["W"])

parameters.set_script("Input",""" 
            # all terms in %....% are replaced by parser
            units metal
            atom_style atomic
            atom_modify map array sort 0 0.0
            read_data  %FirstPathConfiguration% # wildcard
            pair_style hybrid/scaled 0.5 eam/fs 0.5 eam/fs
            pair_coeff * * eam/fs 1 systems/EAM-VAC-W/W.eam.fs W
            pair_coeff * * eam/fs 2 systems/EAM-VAC-W/W.eam.fs W
            run 0
            thermo 10
            run 0
""")

# TEST VALUES
parameters.axes["Temperature"] = [0.,1000.,2000.]
parameters.set("nRepeats",1)

parameters.set("OverDamped",0) # Brownian or Langevin
parameters.set("SampleSteps",2000)
parameters.set("ThermSteps",2000)
parameters.set("ThermWindow",100)

manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)
manager.run()
manager.close()
