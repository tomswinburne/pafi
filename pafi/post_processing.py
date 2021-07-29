"""
raw_ensemble_output has 1 + 4*nWorker colums:

Column 0 : Number of workers with max_jump < user threshold at run time
-> nWorker+1 : av(dF), the gradient
-> 2*nWorker+1 : std(dF), the *worker* variance, *not* ensemble error, should not be zero even in infinite time limit!
-> 3*nWorker+1 : av(Psi) = av(dX).dX_0/|dX_0|^2, the tangent projection
-> 4*nWorker+1 : the maximum per-atom jump following the MD
"""

import io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
import pandas as pd 
import pathlib

kb = 8.617e-5

class Profile():
    """Class for handling (free) energy profiles.
    
    Parameters
    ----------

            data : np.ndarray
                A Nx3 shaped array that contains [reaction coordinate, free energy, and error] as columns.

    Attributes
    ----------

        data : np.ndarray
            The Nx3 array passed as argument.
        barrier : float
            The maximum Free Energy value.
        error : float
            error = error(R=0)+error(R=Rm), where Rm is the reaction corrdinate at which the Free Energy is maximal.
    """

    def __init__(self, data):
        self.data = data
        self.barrier = data[:,1].max()
        self.error = data[0][2] + data[data[:,1].argmax()][2]

class PafiResult():
    """Class for manipulating PAFI simulation results.
    
    Parameters
    ----------

            file : pathlib.Path object or str
                The path to the pafi raw_ensemble file to extract data from. 
                Should have a name of the form "raw_ensemble_output_*K_*", unless optional temperature parameter is given.
            temperature : int or float
                The temperature of the simulation.

    Attributes
    ----------

        temperature : str
            The standard output passed by PreciSo during its execution.
        file : preciso.Distribution object
            A preciso.Distribution object, that contains all the data related to precipitates size distributions.
        discrete_profile : pafi.Profile
            The free energy profile as a pafi.Profile object. Contains the discrete data.
        splined_profile : pafi.Profile
            The free energy profile as a pafi.Profile object. Rediscretized in order to have a smooth profile curve.
        bar :   np.array
            A convenient array that contains commonly used data : [temperature, discrete_profile.barrier, discrete_profile.error, splined_profile.barrier, splined_profile.error]
            See documentation of pafi.Profile for details on .barrier and .error attributes.
    
    Examples
    --------

        >>> import pafi
        >>> import matplotlib.pyplot as plt
        >>> p = pafi.PafiResult('raw_ensemble_output_0K_0')
        >>> ax = p.plot()
        >>> plt.show()
    """

    def __init__(self, file, temperature=None):
        assert (isinstance(file, str)) or (isinstance(file, pathlib.Path)) ,"`file` must be str or pathlib.Path object."

        self.file = file
        if temperature is not None:
            self.temperature = temperature
        else:
            self.temperature = self.get_temperature()

        raw = self.raw_parser()
        self.discrete_profile = Profile(PafiResult.integrate(raw))
 
        remeshed = PafiResult.remesh(raw)
        self.splined_profile = Profile(PafiResult.integrate(remeshed))

        self.bar = self.get_bar()

    def raw_parser(self, disp_thresh=0.7):
        """Parse the raw data file"""
        try:
            f = open(self.file,'r')
        except IOError:
            print("raw data file %s not found" % self.file)
            return np.zeros((1,3))

        count_data = count_valid = 0
        r_dFave_dFstd = []
        r = 0.0
        for line in f.readlines():
            if line[0] != "#":
                fields = np.loadtxt(io.StringIO(line.strip()))
                n_data = (fields.size-1)//4
                n_valid = (fields[-n_data:] < disp_thresh).sum()

                count_data  += n_data
                if n_valid>0:
                    count_valid += n_valid
                
                    r_dFave_dFstd += [[r, 
                                        fields[1:n_valid+1].mean(),
                                        fields[1:n_valid+1].std()/np.sqrt(n_valid)]]
                r += 1.0 # always increment even if n_valid==0
        # print(self.file, " data :", count_data, " valid : ", count_valid, " ({}%)".format(int(count_valid/count_data*100)))
        
        r_dFave_dFstd = np.r_[r_dFave_dFstd]
        r_dFave_dFstd[:,0]/=r_dFave_dFstd[-1][0] # r : 0 -> 1

        return r_dFave_dFstd
    
    def get_bar(self):
        """This data format is used by some legacy functions"""
        return np.array([self.temperature, 
                            self.discrete_profile.barrier, 
                            self.discrete_profile.error, 
                            self.splined_profile.barrier, 
                            self.splined_profile.error])

    def get_temperature(self):
        """The temperature of the run"""
        return int(str(self.file).split("_")[-2][:-1])

    def remesh(data, density=10):
        """Remesh to obtained the splined version of the profile"""
        spl_data = np.zeros((density*data.shape[0],data.shape[1]))
        r_data = np.linspace(0.,1.,data.shape[0])
        r_spl_data = np.linspace(0.,1.,spl_data.shape[0])
        spl_data[:,0] = r_spl_data
        spl_data[:,1] = interp1d(data[:,0], data[:,1],kind='linear',fill_value="extrapolate")(spl_data[:,0])
        spl_data[:,2] = interp1d(data[:,0], data[:,2],kind='linear',fill_value="extrapolate")(spl_data[:,0])
        return spl_data

    def integrate(data, discretization_error_estimate=0.015):
        idata = np.zeros(data.shape)
        idata[:,0] = data[:,0]
        idata[1:,1] = -cumtrapz(data[:,1],data[:,0])
        idata[1:,2] = cumtrapz(data[:,2],data[:,0]) +  discretization_error_estimate
        run_min =  np.minimum.accumulate(idata[:,1])
        run_min_shift =  np.minimum.accumulate(np.append(idata[1:,1],idata[-1][1]))
        if (run_min == run_min_shift).sum()>0:
            idata = idata[run_min == run_min_shift,:]
        idata[:,1]-=idata[0][1]
        return idata

    def plot(self, ax=None, color=None):
        
        if ax is None:
            fig, ax = plt.subplots()
        if color is None:
            color="C0"
        discrete = self.discrete_profile.data
        splined = self.splined_profile.data

        for i, profile in enumerate([self.discrete_profile.data, self.splined_profile.data]):
            
            style = ["o", "-"][i]
            label = [f'{self.temperature}K', ''][i]
            ax.plot(profile[:,0],
                        profile[:,1],
                        style, 
                        color=color,
                        label=label)

            ax.fill_between(profile[:,0],
                    profile[:,1]-profile[:,2],
                    profile[:,1]+profile[:,2],
                    facecolor='0.8')

        ax.set_xlabel("Reaction Coordinate")
        ax.set_ylabel("Free energy of activation [eV]")
        ax.set_xlabel("Temperature [K]")      
        ax.legend(loc="best")  
        return ax

def free_energy_vs_temperature(flist, harmonic_until_T=100, ymin=-0.05):

    fig, ax = plt.subplots(1,2, figsize=(8,4),dpi=144,sharey=True)
    
    data = [PafiResult(file) for file in flist]
    data.sort(key=lambda x: x.temperature)
    barriers = [x.bar for x in data]

    for i, r in enumerate(data):
        a = r.plot(ax[0], color=f'C{i}')

    bar = np.array(barriers)
    bar = bar[np.argsort(bar[:, 0])]

    # Fit a linear regression on the domain below harmonic_until_T K
    harmonic_regime = bar[np.where((bar[:,0]<=harmonic_until_T))]
    p = np.polyfit(harmonic_regime[:,0], harmonic_regime[:,1],1)

    ax[1].plot(bar[:,0], bar[:,1], "-o")
    ax[1].fill_between(bar[:,0],bar[:,1]-bar[:,2],bar[:,1]+bar[:,2],facecolor='0.93')
    ax[1].fill_between(bar[:,0],bar[:,3]-bar[:,4],bar[:,3]+bar[:,4],facecolor='0.93')
    # harmonic domain
    ax[1].plot(bar[:4,0],p[1]+p[0]*bar[:4,0],'k--', label=r'$\Delta U_0=%2.3geV,\Delta S_0=%2.3g{\rm k_B}$' % (p[1], -p[0]/kb))
    ax[1].set_xlabel("Temperature [K]")
    ax[1].legend(loc='best')
    ax[1].set_ylim(ymin=ymin)
    plt.subplots_adjust(wspace=0)
    plt.tight_layout()
        
    return ax