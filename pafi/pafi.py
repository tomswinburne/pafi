"""
raw_ensemble_output has i + 4*nWorker colums. i = 2 since commit 3ca8f884 (20/09/21). For older versions of PAFI, i=1.

Column 0 : r
Column 1 : Number of workers with max_jump < user threshold at run time
-> nWorker+i : av(dF), the gradient
-> 2*nWorker+i : std(dF), the *worker* variance, *not* ensemble error, should not be zero even in infinite time limit!
-> 3*nWorker+i : av(Psi) = av(dX).dX_0/|dX_0|^2, the tangent projection
-> 4*nWorker+i : the maximum per-atom jump following the MD
"""

import io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
import pandas as pd 
import pathlib
from scipy.interpolate import UnivariateSpline

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

    def __init__(self, data, neb_rcoord=None):
        self.data = data
        self.barrier = data[:,1].max() - data[:data[:,1].argmax()].min()
        self.error = data[0][2] + data[data[:,1].argmax()][2]
        self.neb_rcoord = neb_rcoord

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
        file : raw_results file
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

    def __init__(self, file, logfile=None, temperature=None, neb_csv=None, full_path=False):
        assert (isinstance(file, str)) or (isinstance(file, pathlib.Path)) ,"`file` must be str or pathlib.Path object."
        
        self.file = file
        if neb_csv is not None:
            self.neb_rcoord = self.get_r_from_neb_csv(neb_csv)
        else:
            self.neb_rcoord = None
        
        self.r_coord = None
        if logfile is not None:
            self.r_coord = self.log_parser(logfile)
        
        if temperature is not None:
            self.temperature = temperature
        else:
            self.temperature = self.get_temperature()

        raw, self.hist_data, self.std = self.raw_parser()

        self.discrete_profile = Profile(PafiResult.integrate(raw, full_path=full_path), self.neb_rcoord)
 
        remeshed = PafiResult.remesh(raw)
        self.splined_profile = Profile(PafiResult.integrate(remeshed, full_path=full_path), self.neb_rcoord)

        self.bar = self.get_bar()

    def get_r_from_neb_csv(self, file):
        # a = np.genfromtxt(file)
        # a = a[np.isfinite(a).any(axis=1)][:,0]
        # return (a-a[0])/a[-1]

        return np.loadtxt(file)[:,0]
        
    def raw_parser(self, disp_thresh=1.0, mask_outliers=False):
        """Parse the raw data file"""
        try:
            f = open(self.file,'r')
        except IOError:
            print("raw data file %s not found" % self.file)
            return np.zeros((1,3))

        count_data = count_valid = 0
        r_dFave_dFstd = []
        dF_raw = []
        std_raw = []
        

        if PafiResult.has_old_data_format(self.file):
            if self.r_coord is not None:
                # Use rcoord from logfile
                print('Use rcoord from logfile')
                for line, r_c in zip(f.readlines(), self.r_coord):
                    if line[0] != "#":
                        fields = np.loadtxt(io.StringIO(line.strip()))
                        n_data = (fields.size-1)//4
                        n_valid = (fields[-n_data:] < disp_thresh).sum()

                        count_data  += n_data
                        if n_valid>0:
                            count_valid += n_valid
                        
                            r_dFave_dFstd += [[r_c, 
                                                fields[1:n_valid+1].mean(),
                                                fields[1:n_valid+1].std()/np.sqrt(n_valid)]]
                # print(self.file, " data :", count_data, " valid : ", count_valid, " ({}%)".format(int(count_valid/count_data*100)))
                r_dFave_dFstd = np.r_[r_dFave_dFstd]
                r_dFave_dFstd[:,0]/=r_dFave_dFstd[-1][0] # r : 0 -> 1
                return r_dFave_dFstd, len(r_dFave_dFstd)
            else:
                # reaction coordinate is interpolated (may cause deviations at large T)
                # needed for backward compatibility
                print('reaction coordinate is interpolated (may cause deviations at high T)')
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
                return r_dFave_dFstd, len(r_dFave_dFstd)
        else:
            # reaction coordinate is read from file
            print("Reaction coordinate is read from raw file")
            outliers = 0.5
            for line in f.readlines():
                if line[0] != "#":
                    fields = np.loadtxt(io.StringIO(line.strip()))
                    n_data = (fields.size-2)//4
                    n_valid = (fields[-n_data:] < disp_thresh).sum()

                    count_data  += n_data
                    if n_valid>0:
                        count_valid += n_valid

                        dF = np.array(fields[2:n_valid+2])
                        low = np.outer(np.percentile(dF, outliers, axis=0),np.ones(dF.shape[0]))
                        high = np.outer(np.percentile(dF, 100-outliers, axis=0),np.ones(dF.shape[0]))
                        mask = (dF>low) * (dF<high) 
                        masked = (dF*mask).reshape(dF.shape[0])
                        if np.allclose(masked, 0):
                            masked = dF # disable mask for 0K runs
                        if not mask_outliers:
                            masked = dF
                        
                        n = np.count_nonzero(masked)
                        r_dFave_dFstd += [[fields[0],
                                           masked.mean(),
                                           masked.std()/np.sqrt(n)]]

                        dF_raw += list((masked-masked.mean())/masked.std())
                        std_raw += [masked.std()/np.sqrt(n)]
                        # dF_raw += list((fields[2:n_valid+2]-fields[2:n_valid+2].mean())/fields[2:n_valid+2].std())
                        # std_raw += [fields[2:n_valid+2].std()/np.sqrt(n_valid)]
                        # 
                        # r_dFave_dFstd += [[fields[0],
                                         # fields[2:n_valid+2].mean(),
                                         # fields[2:n_valid+2].std()/np.sqrt(n_valid)]]
            
            print(self.file, " data :", count_data, " valid : ", count_valid, " ({}%)".format(int(count_valid/count_data*100)))
            r_dFave_dFstd = np.r_[r_dFave_dFstd]

            return r_dFave_dFstd, dF_raw, std_raw

    def log_parser(self, logfile):
        import re
        with open(logfile, "r") as f:
            content = f.read()
            res = re.findall("""##.+\n[\s\t]*([\d\.]+)""", content)
            R_coord = np.r_[[float(r) for r in res]]
        return R_coord

    def has_old_data_format(file):
        d = np.loadtxt(file, max_rows=1)
        if (d.shape[0]-1) //4 == (d.shape[0]-1) /4:
            return True
        elif (d.shape[0]-2) //4 == (d.shape[0]-2) /4:
            return False
        else: 
            raise RuntimeError(f"File {file} has an unsupported format")

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

    def remesh(data, density=6):
        """Remesh to obtain the splined version of the profile"""
        spl_data = np.zeros((density*data.shape[0], data.shape[1]))
        r_data = np.linspace(0., 1., data.shape[0])
        r_spl_data = np.linspace(0., 1., spl_data.shape[0])
        spl_data[:,0] = interp1d(r_data, data[:,0],kind='linear')(r_spl_data)
        spl_data[:,1] = interp1d(data[:,0], data[:,1],kind='linear',fill_value="extrapolate")(spl_data[:,0])
        spl_data[:,2] = interp1d(data[:,0], data[:,2],kind='linear',fill_value="extrapolate")(spl_data[:,0])
        return spl_data

    def integrate(data, full_path, discretization_error_estimate=1e-4):
        idata = np.zeros(data.shape)
        idata[:,0] = data[:,0]
        idata[1:,1] = -cumtrapz(data[:,1],data[:,0])
        idata[1:,2] = cumtrapz(data[:,2],data[:,0]) +  discretization_error_estimate
        if not full_path:
            run_min =  np.minimum.accumulate(idata[:,1])
            run_min_shift =  np.minimum.accumulate(np.append(idata[1:,1],idata[-1][1]))
            if (run_min == run_min_shift).sum()>0:
                idata = idata[run_min == run_min_shift,:]
        idata[:,1]-=idata[0][1]
        return idata

    def plot(self, ax=None, color=None, label_prepend="", use_neb_rcoord=False):
        
        if ax is None:
            fig, ax = plt.subplots()
        if color is None:
            color="C0"
        if use_neb_rcoord:
            if self.neb_rcoord is None: 
                raise ValueError("Provide neb_csv file with NEB r coordinate as 1st column") 
        discrete = self.discrete_profile.data
        splined = self.splined_profile.data

        for i, profile in enumerate([self.discrete_profile.data, self.splined_profile.data]):
            
            style = ["o", "-"][i]
            label = [label_prepend + f'{self.temperature}K', ''][i]

            if use_neb_rcoord:
                x = np.linspace(0, 1, self.neb_rcoord.shape[0])
                spl = UnivariateSpline(x, self.neb_rcoord, s=0)
                r = spl(np.linspace(0,1, profile[:,1].shape[0]))
                print("using rcorrdinate from NEB file")
            else:
                r = profile[:,0]

            ax.plot(r,
                    profile[:,1],
                    style, 
                    color=color,
                    label=label)
            if i==1: # only errorbars of splined profile
                ax.fill_between(r,
                        profile[:,1]-profile[:,2],
                        profile[:,1]+profile[:,2],
                        facecolor='0.8')

        ax.set_xlabel("Reaction Coordinate")
        ax.set_ylabel("Gibbs energy profile (eV)")
        ax.legend(loc="best")  
        return ax

def free_energy_vs_temperature(flist, ax=None, fit_harmonic=True, harmonic_until_T=100, ymin=-0.05, label_prepend="", 
                                start_color_at_index=0, return_poly=False, add_pts=None, full_path=False):

    if (ax is None) or len(ax)!=2:
        fig, ax = plt.subplots(1,2, figsize=(9,5),dpi=144,sharey=True)
    
    data = [PafiResult(file, full_path=full_path) for file in flist]
    data.sort(key=lambda x: x.temperature)
    bar = np.array([x.bar for x in data])

    # to add arbitrary pts to the right panel plot 
    if add_pts:
        pts = [[pt[0], pt[1], 1e-4, pt[1], 1e-4] for pt in add_pts]
        bar = np.concatenate((np.array(pts), bar))
    for i, r in enumerate(data):
        a = r.plot(ax[0], color=f'C{i+start_color_at_index}', label_prepend=label_prepend)

    # bar = bar[np.argsort(bar[:, 0])]
    bar = bar.astype(float)
    
    if fit_harmonic:
        # Fit a linear regression on the domain below harmonic_until_T K
        harmonic_regime = bar[np.where((bar[:,0]<=harmonic_until_T))]
        p = np.polyfit(harmonic_regime[:,0], harmonic_regime[:,1], 1)

    ax[1].plot(bar[:,0], bar[:,3], "-o", color=f'C{start_color_at_index}', label=label_prepend)
    ax[1].fill_between(bar[:,0], bar[:,1]-bar[:,2], bar[:,1]+bar[:,2],facecolor='0.93')
    ax[1].fill_between(bar[:,0], bar[:,3]-bar[:,4], bar[:,3]+bar[:,4],facecolor='0.93')
    if fit_harmonic:
        # harmonic domain
        ax[1].plot(bar[:,0],p[1]+p[0]*bar[:,0],'--', color=f'C{start_color_at_index}', label=label_prepend + r'$\Delta U_0=%2.3geV,\Delta S_0=%2.3g{\rm k_B}$' % (p[1], -p[0]/kb))
    ax[1].set_xlabel("Temperature (K)")
    ax[1].set_ylabel("Gibbs energy of activation (eV)")
    ax[1].legend(loc='best')
    # ax[1].set_ylim(ymin=ymin)
    plt.subplots_adjust(wspace=0)
    plt.tight_layout()
    
    if return_poly: return ax, p
    else: return ax
