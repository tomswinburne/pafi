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


discretization_error_estimate = 0.005 # estimate in eV, divided by 2

# wildcard for raw dump files
fl = glob.glob("dumps/raw*");


# read in file line-by-line, returning ave and std for data with max displacement <= disp_thresh
def raw_parser(file_name,disp_thresh = 0.4):
    try:
        f = open(file_name,'r')
    except IOError:
        print("raw data file %s not found" % file_name)
        return np.zeros((1,3))

    r_dFave_dFstd = []
    for line in f.readlines():
        if line[0] != "#":
            fields = np.loadtxt(io.StringIO(line.strip()))
            n_data = (fields.size-2)//4
            n_valid = (fields[-n_data:] < disp_thresh).sum()
            if n_valid>0:
                r_dFave_dFstd += [[fields[0],fields[2:n_valid+2].mean(),fields[2:n_valid+2].std()/np.sqrt(n_valid)]]
    r_dFave_dFstd = np.r_[r_dFave_dFstd]
    # r_dFave_dFstd[:,0]/=r_dFave_dFstd[-1][0] # r : 0 -> 1

    return r_dFave_dFstd

# splined rediscretization
def remesh(data,density = 10):
    spl_data = np.zeros((density*data.shape[0],data.shape[1]))
    r_data = np.linspace(0.,1.,data.shape[0])
    r_spl_data = np.linspace(0.,1.,spl_data.shape[0])
    spl_data[:,0] = interp1d(r_data, data[:,0],kind='linear')(r_spl_data)
    spl_data[:,1] = interp1d(data[:,0], data[:,1],kind='linear')(spl_data[:,0])
    spl_data[:,2] = interp1d(data[:,0], data[:,2],kind='linear')(spl_data[:,0])
    return spl_data

def integrate(data,find_first_min=False):
    idata = np.zeros(data.shape)
    idata[:,0] = data[:,0]
    idata[1:,1] = -cumtrapz(data[:,1],data[:,0])
    idata[1:,2] = cumtrapz(data[:,2],data[:,0]) +  discretization_error_estimate
    if find_first_min:
        run_min =  np.minimum.accumulate(idata[:,1])
        run_min_shift =  np.minimum.accumulate(np.append(idata[1:,1],idata[-1][1]))
        if (run_min == run_min_shift).sum()>0:
            idata = idata[run_min == run_min_shift,:]
    idata[:,1]-=idata[0][1]
    return idata, idata[:,1].max(), idata[0][2] + idata[idata[:,1].argmax()][2]



# temperatures
T = np.r_[[int(f.split("_")[-2][:-1]) for f in fl]];


fig,axs = plt.subplots(1,2,figsize=(8,4),dpi=144,sharey=True);

bar = []
for ii,i_f in enumerate(T.argsort()):
    _bar = [T[i_f]]
    r_dFave_dFstd = raw_parser(fl[i_f])
    r_Fave_Fstd,barrier,error = integrate(r_dFave_dFstd)

    _bar += [barrier]+[error]
    r_dFave_dFstd_remesh = remesh(r_dFave_dFstd)
    r_Fave_Fstd_remesh,barrier,error = integrate(r_dFave_dFstd_remesh)
    _bar += [barrier]+[error]

    axs[0].plot(r_Fave_Fstd[:,0],
                r_Fave_Fstd[:,1],
                f'C{ii%9}o--',label='%dK' % T[i_f])

    axs[0].fill_between(r_Fave_Fstd[:,0],
            r_Fave_Fstd[:,1]-r_Fave_Fstd[:,2],
            r_Fave_Fstd[:,1]+r_Fave_Fstd[:,2],
            facecolor='0.8')

    axs[0].plot(r_Fave_Fstd_remesh[:,0],
                r_Fave_Fstd_remesh[:,1],
                f'C{ii%9}-',label='%dK (Splined)' % T[i_f])

    axs[0].fill_between(r_Fave_Fstd_remesh[:,0],
            r_Fave_Fstd_remesh[:,1]-r_Fave_Fstd_remesh[:,2],
            r_Fave_Fstd_remesh[:,1]+r_Fave_Fstd_remesh[:,2],
            facecolor='0.8')

    bar += [_bar]
bar = np.r_[bar]
print(bar.shape)

p = np.polyfit(bar[:3,0],bar[:3,1],1)

kb = 8.617e-5

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
