from mpi4py import MPI
import os
import numpy as np
from typing import List
from ..parsers.PAFIParser import PAFIParser
from .ResultsHolder import ResultsHolder

class BaseGatherer:
    def __init__(self,params:PAFIParser,
                 nWorkers:int,
                 rank:int,
                 ensemble_comm:MPI.Intracomm,
                 roots:List[int])->None:
        """Basic gatherer of PAFI simulation data

        Parameters
        ----------
        params : PAFIParser
            Custom or predefined PAFIParser object
        nWorkers : int
            total number of PAFI workers
        rank : int
            global MPI rank
        ensemble_comm : MPI.Intracomm
            MPI communicator to gather ensemble data
        roots : List[int]
            list of root ranks for each worker. Only collate data here
            This should be depreceated asap
        """
        self.params = params
        self.nWorkers = nWorkers
        self.rank = rank
        self.comm = ensemble_comm
        self.roots = roots
        self.epoch_data = None # for each cycle
        self.all_data = None # for total simulation
        self.last_data = None # for print out
    
    def gather(self,data:dict|ResultsHolder)->None:
        """Gather results from a simulation epoch,
        local to each worker. Here, very simple,
        is overwritten by each call of gather()

        Parameters
        ----------
        data : dict or ResultsHolder
            Simulation data, extracted as dictionary from ResultsHolder
        """
        if self.rank in self.roots:
            if isinstance(data,ResultsHolder):
                self.epoch_data = data.data.copy()
            else:
                self.epoch_data = data.copy()
       
    
    def collate(self,repeat:int=0)->None:
        """Collate all data on root node
        It is assumed that all data is in dictionary form, with identical keys
        Only basic checks for multiple calls.

        Parameters
        ----------
        repeat : int, optional
            If  0, wipe last_data. Default 0
        """
        if self.rank in self.roots:
            if not self.epoch_data is None:
                all_epoch_data = self.comm.gather(self.epoch_data)
                self.epoch_data = None
            else:
                all_epoch_data = None

            if self.rank == 0 and not all_epoch_data is None:
                if (not repeat) or self.last_data is None:
                    self.last_data = {k:[] for k in all_epoch_data[0].keys()}
                if repeat:
                    self.last_data["Repeat"] = []

                if self.all_data is None:
                    self.all_data = {k:[] for k in all_epoch_data[0].keys()}
                
                for k in all_epoch_data[0].keys():
                    self.all_data[k] += list(d[k] for d in all_epoch_data)
                    self.last_data[k] += list(d[k] for d in all_epoch_data)

    
    def get_line(self,fields:List[str])->List[str]|None:
        """Return output data to print out on root node

        Parameters
        ----------
        fields : List[str]
            fields to extract
        Returns
        -------
        List[str]|None
            if root process, return list of lines to print, else return `None`
        """
        if self.rank != 0:
            return None
        line = []
        for f in fields:
            std = bool(f[-4:] == "_std")
            key = f[:-4] if std else f
            if (not self.last_data is None) and (key in self.last_data.keys()):
                d = self.last_data[key]
                line += [np.std(d)/np.sqrt(len(d)) if std else np.mean(d)]
            else:
                line += ["n/a"]
            
        return line
    
    def get_dict(self,fields:List[str])->dict|None:
        """Return output data to print out on root node

        Parameters
        ----------
        fields : List[str]
            fields to extract
        Returns
        -------
        dict|None
            if root process, return dict of fields to print, else return `None`
        """
        line = self.get_line(fields)
        if not line is None:
            return {kv[0]:kv[1] for kv in zip(fields,line)}
    
    def write_pandas(self,path:os.PathLike[str])->None:
        """Write data as pandas dataframe

        Parameters
        ----------
        path : os.PathLike[str]
            path to file
        """

        if self.rank==0:
            if self.all_data is None:
                print("No data to write! Exiting!")
            else:
                import pandas as pd
                if not os.path.isdir(os.path.dirname(path)):
                    raise IOError("Unknown directory for writing csv!")
                else:
                    df = pd.DataFrame(self.all_data)
                    df.metadata = " ".join(list(self.params.axes.keys()))
                    #df.attrs = self.params.to_dict().copy()
                    df.to_csv(path)
    
    def read_pandas(self,path:os.PathLike[str])->None:
        """Read in data from pandas dataframe

        Only possible if no sampling has taken place.

        Parameters
        ----------
        path : os.PathLike[str]
            path to file
        """
        assert self.all_data is None
        assert os.path.exists(path)
        import pandas as pd
        self.all_data = pd.read_csv(path)
            
                
