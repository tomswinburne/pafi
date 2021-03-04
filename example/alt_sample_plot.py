"""
VERY rough Python code to plot PAFI results

(c) 2021 Tom Swinburne

raw_ensemble_output has 1 + 4*nWorker colums:

Column 0 : Number of workers with max_jump < user threshold at run time
-> nWorker+1 : av(dF), the gradient
-> 2*nWorker+1 : std(dF), the *worker* variance, *not* ensemble error, should not be zero even in infinite time limit!
-> 3*nWorker+1 : av(Psi) = av(dX).dX_0/|dX_0|^2, the tangent projection
-> 4*nWorker+1 : the maximum per-atom jump following the MD
"""

import glob,io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

# read in file line-by-line, returning ave and std for data with max displacement <= disp_thresh
def raw_parser(file_name,disp_thresh = 0.4):
    try:
        f = open(file_name,'r')
    except IOError:
        print("raw data file %s not found" % file_name)
        return np.zeros((1,3))

    r_dFave_dFstd = []
    r = 0.0
    for line in f.readlines():
        if line[0] != "#":
            fields = np.loadtxt(io.StringIO(line.strip()))
            n_data = (fields.size-1)//4
            n_valid = (fields[-n_data:] < disp_thresh).sum()
            if n_valid>0:
                r_dFave_dFstd += [[r,fields[1:n_valid+1].mean(),fields[1:n_valid+1].std()/np.sqrt(n_valid)]]
            r += 1.0 # always increment even if n_valid==0
    r_dFave_dFstd = np.r_[r_dFave_dFstd]
    r_dFave_dFstd[:,0]/=r_dFave_dFstd[-1][0] # r : 0 -> 1

    return r_dFave_dFstd



# file list- here we take epoch 5
fl = glob.glob("data/50/raw_ensemble_output_*_5");

# temperatures
T = np.r_[[int(f.split("_")[-2][:-1]) for f in fl]];

disp_thresh = 0.4 # new max jump threshold

nReSpl = 10 # resplining density

discretization_error_estimate = 0.015 # estimate in eV, divided by 2


fig,axs = plt.subplots(1,2,figsize=(8,4),dpi=144,sharey=True);


bar = []
for ii,i_f in enumerate(T.argsort()):
    rdFdFe = raw_parser(fl[i_f])


    spl_rdFdFe = np.zeros((nReSpl*rdFdFe.shape[0],rdFdFe.shape[1]))
    spl_rdFdFe[:,0] = np.linspace(0.,1.,spl_rdFdFe.shape[0])
    spl_rdFdFe[:,1] = interp1d(rdFdFe[:,0], rdFdFe[:,1],kind='linear')(spl_rdFdFe[:,0])
    spl_rdFdFe[:,2] = interp1d(rdFdFe[:,0], rdFdFe[:,2],kind='linear')(spl_rdFdFe[:,0])

    _bar = [T[i_f]]
    for j,_rdFdFe in enumerate([rdFdFe,spl_rdFdFe]):
        rFFe = np.zeros(_rdFdFe.shape)
        rFFe[:,0] = np.linspace(0.,1.,_rdFdFe.shape[0])
        rFFe[1:,1] = -cumtrapz(_rdFdFe[:,1],_rdFdFe[:,0])
        rFFe[1:,2] = cumtrapz(_rdFdFe[:,2],_rdFdFe[:,0]) +  discretization_error_estimate

        # integrate from minimum (important for barrier under stress)
        first_argmin = 0
        while rFFe[first_argmin][1] == rFFe[:first_argmin+1,1].min():
            first_argmin+=1
        first_argmin-=1
        rFFe = rFFe[first_argmin:,:]

        rFFe[:,0] -= rFFe[0][0]
        rFFe[:,1] -= rFFe[:,1][0]




        axs[0].plot(rFFe[:,0],rFFe[:,1],'C%d%s' % (ii,['o--','-'][j]),label='%dK%s' % (T[i_f],[""," (Splined)"][j]))
        axs[0].fill_between(rFFe[:,0],rFFe[:,1]-rFFe[:,2],rFFe[:,1]+rFFe[:,2],facecolor='0.8')

        barrier = rFFe[:,1].max()-rFFe[:,1].min()
        error_barrier = rFFe[:,2][rFFe[:,1].argmax()]+rFFe[:,2][rFFe[:,1].argmin()]
        _bar += [barrier]+[error_barrier]
    bar += [_bar]
bar = np.r_[bar]

p = np.polyfit(bar[:3,0],bar[:3,1],1)

kb = 8.617e-5
b = 2.8552 * np.sqrt(.75)

axs[1].plot(bar[:,0],bar[:,1],'o-',label='Raw Data')
axs[1].plot(bar[:,0],bar[:,3],'o--',label='Splined Data')

axs[1].plot(bar[:,0],p[1]+p[0]*bar[:,0],'k--',label=r'$\Delta U_0=%2.3geV,\Delta S_0=%2.3g{\rm k_B}$' % (p[1],-p[0]/kb))
axs[1].fill_between(bar[:,0],bar[:,1]-bar[:,2],bar[:,1]+bar[:,2],facecolor='0.8')
axs[1].fill_between(bar[:,0],bar[:,3]-bar[:,4],bar[:,3]+bar[:,4],facecolor='0.8')



axs[1].set_ylim(ymin=-discretization_error_estimate)
axs[0].legend(fontsize=7)
axs[1].legend(fontsize=7)

axs[0].set_xlabel("Reaction Coordinate")
axs[1].set_xlabel("Temperature [K]")

axs[0].set_ylabel("Free energy barrier [eV]")

plt.subplots_adjust(wspace=0)
plt.savefig("PAFI_fig.pdf")
