from mpi4py import MPI
from pafi import PAFIManager

# Read in all parameters from configuration file
config = "./configuration_files/CompleteConfiguration_TEST.xml"
# or overwrite default parameters with partial configuration file
# config = "./configuration_files/PartialConfiguration_TEST.xml"

manager = PAFIManager(MPI.COMM_WORLD,config)
manager.run()
manager.close()