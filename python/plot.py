import pafi
import sys
import matplotlib.pyplot as plt
from glob import glob
import numpy as np


"""
Rough plotting script that calls pafi.py from the command line
Todo: port to click CLI / argparse / other

Usage 

Plot a single file (pass "s" as 1st arg)
python3 plot.py s dumps/raw_ensemble_output_40K_0

Plot a single file, projecting on neb reaction coordinate
python3 plot.py s dumps/raw_ensemble_output_40K_0 rneb neb_coord.csv

Plot a series a runs using a glob() expression (pass "m" as 1st arg)
python3 plot.py m "T*/pafi/raw*"

Add a point at (0, 0.4343) obtained from NEB (keyword "zero")
python3 plot.py m "T*/pafi/raw*" zero 0.4343
"""

print(sys.argv)
    
fig,axs = plt.subplots(1,3, figsize=(9,3))

if sys.argv[-2] == "rneb":
    # to take r from neb if RealMEP=0
    neb_coord = sys.argv[-1]
else:
    neb_coord = None
    
if sys.argv[1] == "s":
    
    plt.sca(axs[0])
    for file in sys.argv[2:]:

        if file == "rneb": break # hack to deal with command line args

        # plot barriers
        r = pafi.PafiResult(file, neb_csv=neb_coord)
        r.plot(ax=axs[0], use_neb_rcoord=(neb_coord!=None))
    try:
        ## will fail for 0K runs
        plt.sca(axs[1])
        plt.hist(r.hist_data, bins=30)
        plt.ylabel("(dF-dFmean)/std")

        plt.sca(axs[2])
        plt.plot(r.std, "o-")
        plt.ylabel("std")
        plt.xlabel("hyperplane")
    except: pass
    finally:
        plt.tight_layout()
        plt.show()

elif sys.argv[1] == "m":
    # for multiple files

    if len(sys.argv)>=4:
        if sys.argv[3] == "zero":
            pafi.free_energy_vs_temperature(glob(sys.argv[2]), add_pts=[(0, sys.argv[4])])
        elif sys.argv[3] == "neb":
            d = np.loadtxt("neb/barrier_free_end.csv", skiprows=1)
            zero = d[:,1].max()
            print(zero)
            pafi.free_energy_vs_temperature(glob(sys.argv[2]), add_pts=[(0, zero)])
    else:
        pafi.free_energy_vs_temperature(glob(sys.argv[2]))
plt.show()

     
