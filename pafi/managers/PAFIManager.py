import itertools
from typing import List,Dict
import numpy as np
import os
from mpi4py import MPI
from ..results.ResultsHolder import ResultsHolder
from .BaseManager import BaseManager
from ..parsers.PAFIParser import PAFIParser
from ..workers.PAFIWorker import PAFIWorker
from ..results.Gatherer import Gatherer

class PAFIManager(BaseManager):
    def __init__(self, world: MPI.Intracomm, 
                 xml_path:None|os.PathLike[str]=None,
                 parameters:None|PAFIParser=None,
                 restart_data:None|os.PathLike[str]=None,
                 Worker:PAFIWorker=PAFIWorker,
                 Gatherer:Gatherer=Gatherer) -> None:
        """Default manager of PAFI, child of BaseManager

        Parameters
        ----------
        world : MPI.Intracomm
            MPI communicator
        xml_path : None or os.PathLike[str], optional
            path to XML configuration file, default None
        parameters : None or PAFIParser object, optional
            preloaded PAFIParser object, default None
        restart_data : None or os.PathLike[str], optional
            path to CSV data file. Will read and skip already sampled parameters
        Worker : PAFIWorker, optional,
            Can be overwritten by child class, by default PAFIWorker
        Gatherer : Gatherer, optional
            Can be overwritten by child class, by default Gatherer
        """
        
        
        assert (not parameters is None) or (not xml_path is None)
        
        if parameters is None:
            # TODO have standalone check for suffix in config_[suffix].xml?
            # not the best solution currently if we use BaseManager alone....
            parameters = PAFIParser(xml_path=xml_path,rank=world.Get_rank())
        
        super().__init__(world, parameters, Worker, Gatherer)
        
    
    
    def run(self,print_fields:List[str]|None=None,
            width:int=10,
            precision:int=5)->None:
        """Basic parallel PAFI sampling

            Performs a nested loop over all <Axes>, in the order
            presented in the XML configuration file.
            Here, parallelization is naive- all workers are given the 
            same parameters.
        Parameters
        ----------
        print_fields : List[str] or None
            Fields to print to screen, default None. 
            If None, will print "Temperature","ReactionCoordinate","FreeEnergyGradient"
        width : int
            character count of field printout, default 10
        precision : int
            precision of field printout, default 4
        """
        assert self.parameters.ready()

        if print_fields is None:
            print_fields = \
                ["Temperature","postTemperature","ReactionCoordinate","FreeEnergyGradient","FreeEnergyGradient_std"]
        
        nRepeats = 1
        if not self.parameters("nRepeats") is None:
            if self.parameters("nRepeats")>1:
                nRepeats = self.parameters("nRepeats")
                print_fields = ["Repeat"] + print_fields
        
        for f in print_fields:
            width = max(width,len(f))
        
        def line(data:List[float|int|str]|Dict[str,float|int|str])->str:
            """Format list of results to print to screen
            """
            if len(data) == 0:
                return ""
            format_string = ("{: >%d} "%width)*len(data)
            if isinstance(data,dict):
                _fields = []
                for f in print_fields:
                    if f=='Repeat':
                        val = f"{int(data[f])}/{nRepeats}"
                    else:
                        val = data[f]
                    _fields += [val]
            else: 
                _fields = data
            
            fields = []
            for f in _fields:
                isstr = isinstance(f,str)
                fields += [f if isstr else np.round(f,precision)]
            return format_string.format(*fields)

        
        if self.rank==0:
            screen_out = f"""
            Initialized {self.nWorkers} workers with {self.CoresPerWorker} cores
            <> == time averages,  av/err over ensemble
            """
            if min(self.parameters.axes["Temperature"]) < 0.1:
                screen_out+="""
            *** FOR T=0K RUNS SampleSteps=1 AND ThermalSteps=1 ***
            """
            print(screen_out)
            print(line(print_fields))
            
            # return value
            average_results = {k:[] for k in print_fields}
            


        
        last_coord = None
        for axes_coord in itertools.product(*self.parameters.axes.values()):
            dict_axes = dict(zip(self.parameters.axes.keys(), axes_coord))
            
            if self.rank==0:
                if not last_coord is None and last_coord!=axes_coord[:-1]:
                    print("\n"+line(print_fields))
            if nRepeats>1:
                dict_axes["Repeat"] = 1
            results = ResultsHolder()
            results.set_dict(dict_axes)
            
            
            # Useful helper for including zero temperature cheaply...
            for k in ["SampleSteps","ThermSteps","ThermWindow"]:
                if results("Temperature")<0.1:
                    results.set(k,1)
                else:
                    results.set(k,self.parameters(k))
            
            for repeat in range(nRepeats):
                # Sampling run, returning ResultsHolder object
                final_results = self.Worker.sample(results)
                final_results.set("Repeat",repeat + 1)

                # incorporate results (this is only performed on local roots)
                self.Gatherer.gather(final_results)
                self.Gatherer.collate(repeat)

                # wait
                self.world.Barrier()
                screen_out = self.Gatherer.get_dict(print_fields)
                if self.rank == 0:
                    for k in screen_out.keys():
                        average_results[k].append(screen_out[k])
                    
                    print(line(screen_out))
                    self.Gatherer.write_pandas(path=self.parameters.csv_file)
            
            last_coord = axes_coord[:-1]
        
        if self.rank==0:
            print(f"Data written to {self.parameters.csv_file}")
            return average_results
        else:
            return None

        
    