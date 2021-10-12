from pafi.post_processing import *

def cmd_plot():
    import sys
    import matplotlib.pyplot as plt
    p = PafiResult(sys.argv[1])
    ax = p.plot()
    if len(sys.argv)==3:
        plt.savefig(sys.argv[2])
    else:
        plt.show()