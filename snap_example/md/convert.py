from lammps import lammps
import numpy as np
import glob
s = 0.9920217982124345

load = lambda fn :"""units metal
atom_style atomic
atom_modify map array sort 0 0.0
neigh_modify every 2 delay 10 check yes page 1000000 one 100000
read_data  {0}
pair_style    snap
pair_coeff * * Fe.snapcoeff Fe.snapparam Fe
change_box all x scale {1} y scale {1} z scale {1} remap
write_data {2}
""".format(fn,s,"snap_"+fn)

for ff in glob.glob("image_*"):
  lmp=lammps()
  lmp.commands_string(load(ff))
  lmp.close()
