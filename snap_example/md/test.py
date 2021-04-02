from mpi4py import MPI
from lammps import lammps
import numpy as np
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

lmp=lammps()
lmp.file("test.lammps")

#print(f.mean(axis=0))
c = np.ctypeslib.as_array(lmp.gather("c_b",1,30)).reshape((-1,30)).mean(axis=0)
if rank==0:
    f = []
    for i in range(30):
        f += [lmp.extract_fix("ab",0,1,i)]
    print(np.r_[f])
    print("---")
    print(c-np.r_[f])

lmp.close()
MPI.Finalize()
