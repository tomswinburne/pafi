from mpi4py import MPI
from pafi import PAFIManager

complete = True

if complete:
    """
        Load in complete configuration file
    """
    config = "./configuration_files/CompleteConfiguration_TEST.xml"
else:
    """
        Overwrite default parameters with partial configuration file
    """
    config = "./configuration_files/PartialConfiguration_TEST.xml"


manager = PAFIManager(MPI.COMM_WORLD,config)
manager.run()
manager.close()