import pandas as pd
import os,itertools
import numpy as np
from typing import List
from ..parsers.PAFIParser import PAFIParser
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d

class ResultsProcessor:
    def __init__(self,
                 data_path:os.PathLike[str]|List[os.PathLike[str]],
                 xml_path:None|os.PathLike[str]=None,
                 axes:List[str]=None) -> None:
        """Read in PAFI data and plot results

        Parameters
        ----------
        data_path : os.PathLike[str] | List[os.PathLike[str]]
            path, wildcard, or list of paths to PAFI csv files
        xml_path : None | os.PathLike[str]
            path to PAFI XML configuration file, default None
        axes : List[str], optional
            List of axes, overwritten by data in xml_path if present,
            default : None, will be set to ["ReactionCoordinate","Temperature"]

        Methods
        ----------
        append()
        extract_axes()
        integrate()
        """

        self.data = None
        self.axes = None
        self.fields = None
        if not isinstance(data_path,list):
            data_path = [data_path]
        for dp in data_path:
            if self.data is None:
                if os.path.exists(dp):
                    self.data = pd.read_csv(dp)
            else:
                self.append(dp)
        
        if (not xml_path is None) and os.path.exists(xml_path):
            self.params = PAFIParser(xml_path,postprocessing=True)
        else:
            self.params = None
        
        self.extract_axes(axes=axes)
    
    def extract_axes(self,axes=None):
        self.fields = list(self.data.keys())[1:]
        self.axes = {}

        if not self.params is None:
            axes = self.params.axes.keys()
        else:
            if axes is None:
                axes = ["ReactionCoordinate","Temperature"]
            print(f"""
                No config_*.xml file specified (via 'xml_path' argument)
                Using axes (sampling marginals): {axes}
                """)
        for a in axes:
            self.axes[a] = list(np.round(np.sort(np.unique(self.data[a])),4))

    def append(self,data_path)->None:
        if os.path.exists(data_path):
            new_data = pd.read_csv(data_path)
        
        """Append DataFrame to existing stored data. 
        Will throw error if `new_data` does not have fields of `self.data`stored fields to not cover append field

        Parameters
        ----------
        pd : pd.DataFrame
            DataFrame to append
        """
        if set(new_data.keys())>=set(self.data.keys()):
            old_data = self.data.copy()
            self.data = {}
            for f in old_data.keys():
                self.data[f] = \
                    pd.concat([old_data[f],new_data[f]],ignore_index=True)
            self.data = pd.DataFrame(self.data)
        else:
            print("Could not append!")
        
        # recheck ranges
        if not self.axes is None:
            for a in self.axes:
                self.axes[a] = \
                    list(np.round(np.sort(np.unique(self.data[a])),4))

    def ensemble_collate(self,return_pd:bool=False)->None|pd.DataFrame:
        """Perform ensemble averaging 

       
        Each non-axis field F is averaged and returned as 
        F_ave and F_std, using only valid values. 
        Parameters
        ----------
            return_pd : bool, optional
                Return stored pandas dataframe 'ave_data', default False
        Returns
        -------
        pd.DataFrame
            average dataframe, if `return_pd` is True
        """
        valid_key = 'Valid'
        self.count_key = 'ValidCount'
        axes = set(self.axes.keys())
        ave_data = {f:[] for f in self.fields}
        var_data = {f+"_var":[] for f in self.fields}
        count_data = {self.count_key:[]}
        
        # iterate over all axes
        for pt in itertools.product(*[self.axes[a] for a in axes]):
            sel = np.ones_like(self.data[valid_key].to_numpy(),bool)
            for a in zip(axes,pt):
                sel *= np.isclose(self.data[a[0]].to_numpy(),a[1])
            sel_data = self.data[sel]
            
            valid_r = sel_data[valid_key]
            valid_c = valid_r.sum()
            count_data[self.count_key] += [valid_c]

            for f in self.fields:
                if valid_c < 1:
                    ave_data[f] += [0.0]
                    var_data[f+"_var"] += [0.0]
                else:
                    valid_d = sel_data[f][valid_r].to_numpy()
                    ave_data[f] += [valid_d.mean()]
                    var_data[f+"_var"] += [valid_d.var()/valid_c]

        
        self.ave_data = pd.DataFrame({**ave_data,**count_data,**var_data})
        if return_pd:
            return self.ave_data

    def integrate(self,
                  argument:str='ReactionCoordinate',
                  target:str='FreeEnergyGradient',
                  variance:str='FreeEnergyGradientVariance',
                  remesh:int=5,
                  return_remeshed_array:bool=False)->pd.DataFrame:
        
        """Cumulative integration of data along an axis. 
        Integration routine makes a spline interpolation
        to increase the number of intergrand evaluations 
        
        Parameters
        ----------
        argument : str, optional
            integration argument, by default 'ReactionCoordinate'
        target : str, optional
            integral, by default 'FreeEnergyGradient'
        variance : str, optional
            error term, by default 'FreeEnergyGradientVariance'

        remesh : int, optional
            the number of integrand evaluations between 
            existing knot points, by default 5
        return_remeshed_array : bool, optional,
            return dense numpy array for plotting. Default False
        Returns
        -------
            DataFrame with new field '`target`_ave_integrated'
            if return_remeshed_array, also return dense numpy array
        """
        
        # redo ensemble average
        self.ensemble_collate(return_pd=False)
        data = self.ave_data.copy()
        
        # determine keys
        x_key = argument
        y_key = target
        v_key = target+"_var"
        y_key_std = target+"_err"

        i_key = target+"_integrated"
        i_key_std = target+"_integrated_err"
        
        i_keys = [i_key,i_key_std,i_key+"_u",i_key+"_l"]
        i_data = np.zeros((data[x_key].size,4))

        remesh = max(1,remesh)

        # sanity check!
        assert x_key in self.axes
        assert y_key in self.ave_data.keys()
        
        # we will integrate over constant aux values
        auxs = list(set(self.axes) - set([x_key]))
        
        if return_remeshed_array:
            out_array = []
        
        for pt in itertools.product(*[self.axes[a] for a in auxs]):
            # dataframe of selected data
            sel = np.ones_like(data[x_key].to_numpy(),bool)

            for a in zip(auxs,pt):
                sel *= np.isclose(data[a[0]].to_numpy(),a[1])
            sel_data = data[sel]
            

            x_val = np.sort(np.unique(sel_data[x_key]))
            y_val = np.zeros((x_val.size,2))
            for i,x in enumerate(x_val):
                x_sel_data = sel_data[np.isclose(sel_data[x_key],x)]
                y_val[i][0] = x_sel_data[y_key]
                y_val[i][1] = x_sel_data[v_key]

            y_spl = interp1d(x_val,y_val,axis=0,kind='cubic')
            dense_x = np.linspace(x_val.min(),x_val.max(),remesh*x_val.size)
            dense_y = y_spl(dense_x)
            dense_i = cumulative_trapezoid(dense_y,dense_x,axis=0,initial=0)

            i_spl = interp1d(dense_x,dense_i,axis=0,kind='cubic')
            
            if return_remeshed_array:
                out_row = {a[0]:a[1] for a in zip(auxs,pt)}
                out_row[x_key] = dense_x
                out_row[y_key] = dense_y[:,0]
                out_row[v_key] = dense_y[:,1]
                out_row[y_key_std] = np.sqrt(np.abs(dense_y[:,1]))
                out_row[i_key] = dense_i[:,0]
                out_row[i_key_std] = np.sqrt(np.abs(dense_i[:,1]))
                out_array += [out_row]

            for i,x in enumerate(x_val):
                i_dat = i_spl(x)
                i_row = [i_dat[0],i_dat[1],i_dat[0]+i_dat[1],i_dat[0]-i_dat[1]]
                i_data[sel*np.isclose(data[x_key],x)] = i_row
            
            
        for i,k in enumerate(i_keys):
            data[k] = i_data[:,i]
        if return_remeshed_array:
            return data, out_array
        else:
            return data
        """
        # will be converted to dataframe
        out_data = []
        blank_row = {k:[] for k in list(axes)+[y_key,y_key+"_u",y_key+"_l"]}
        if not v_key is None:
            blank_row[v_key] = []
            blank_row[v_key+"_u"] = []
            blank_row[v_key+"_l"] = []
        
        # iterate over all aux axes
        for pt in itertools.product(*[self.axes[a] for a in auxs]):
            
            # select data at constant aux values
            sel = np.ones_like(data[x_key].to_numpy(),bool)
            row = blank_row.copy()
            for a in zip(auxs,pt):
                sel *= np.isclose(data[a[0]].to_numpy(),a[1])
            sel_data = data[sel]
            
                   
            # find integration points, sorted
            x_val = np.sort(np.unique(sel_data[x_key]))
            for a in zip(auxs,pt):
                row[a[0]] = a[1] * np.ones_like(x_val)


            if not v_key is None:
                y_val = np.zeros((x_val.size,4))
            else:
                y_val = np.zeros((x_val.size,2))
            
            for i,x in enumerate(x_val):
                raw_match = sel_data[np.isclose(sel_data[x_key],x)]
                valid_r = raw_match[valid_key]
                y_r = raw_match[y_key][valid_r]
                # see write up for derivation of additional 1/N factor
                y_val[i][:2] = np.r_[[y_r.mean(), y_r.var() / y_r.size]]
                if not v_key is None:
                    v_r = sel_data[np.isclose(sel_data[x_key],x)][v_key]
                    y_val[i][2:] = np.r_[[v_r.mean(), v_r.var() / v_r.size]]

            
            y_spl = interp1d(x_val,y_val,axis=0)
            x_val = np.linspace(x_val.min(),x_val.max(),remesh*x_val.size)
            y_val = y_spl(x_val)
            
            row[x_key] = x_val
            row[y_key] = np.append(0.,cumtrapz(y_val[:,0],x_val))
            row[y_key+"_u"] = np.append(0.,cumtrapz(y_val[:,1],x_val))
            row[y_key+"_u"] = np.sqrt(row[y_key+"_u"])+row[y_key]
            row[y_key+"_l"] = 2.0*row[y_key] - row[y_key+"_u"]
            if not v_key is None:
                row[v_key] = np.append(0.,cumtrapz(y_val[:,2],x_val))
                row[v_key+"_u"] = np.append(0.,cumtrapz(y_val[:,3],x_val))
                row[v_key+"_u"] = np.sqrt(row[v_key+"_u"])+row[v_key]
                row[v_key+"_l"] = 2.0*row[v_key] - row[v_key+"_u"]
            out_data += [row]
        
        return out_data
        """