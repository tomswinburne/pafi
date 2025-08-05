from mpi4py import MPI
from pafi import PAFIManager

# Read in all parameters from configuration file
config = "./xml/CompleteConfiguration_TEST.xml"
# or overwrite default parameters with partial configuration file
# config = "./xml/PartialConfiguration_TEST.xml"

manager = PAFIManager(MPI.COMM_WORLD,config)
manager.run()
manager.close()
