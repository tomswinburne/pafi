import numpy as np
import os
import numpy as np
from typing import Union,Any
ScriptArg = Union[int,float,str]

from .BaseParser import BaseParser

class PAFIParser(BaseParser):
    """Default reader of PAFI XML configuration file
    
    Parameters
    ----------
    xml_path : os.PathLike[str]
        path to XML file, default None
    postprocessing: bool, optional
        are we just looking in postprocess?, default False
    Methods
    ----------
    __call__
    find_suffix_and_write
    

    Raises
    ------
    IOError
        If path is not found
    """
    def __init__(self,
            xml_path:None|os.PathLike[str]=None,
            postprocessing:bool=False,
            rank:int=0) -> None:
        super().__init__(xml_path,postprocessing,rank)

        # initial seed, but must be different across workers...
        self.seeded = False
        
        # repeats, valid bounds (see set_min_valid())
        self.nRepeats = self.parameters["nRepeats"]
        self.maxExtraRepeats = self.parameters["maxExtraRepeats"]         
        self.set_min_valid(1) # temporary
        
    
    def set_min_valid(self,nWorkers:int)->None:
        """Set the minimum number of value samples 
        for each hyperplane

        Parameters
        ----------
        nWorkers : int
            number of MPI workers
        """
        
        self.minValidResults = self.parameters["ReSampleThresh"]
        self.minValidResults *= self.parameters["nRepeats"]
        self.minValidResults *= nWorkers
        self.minValidResults = int(self.minValidResults)
    def set(self,key:str,value:Any,create:bool=False)->None:
        """Set a parameter

        Parameters
        ----------
        key : str
            key of parameter. Must already exist if create is `False`
        value : Any
            value for entry
        create : bool, optional
            Create new entry if true    
        """
        if not create:
            assert key in self.parameters.keys()
        self.parameters[key] = value
    
    def min_valid(self)->int:
        """Return the minimum number of value samples, 
        set via set_min_valid()

        Returns
        -------
        int
            the return value
        """
        return self.minValidResults
    
    def seed(self,worker_instance:int)->None:
        """Generate random number seed

        Parameters
        ----------
        worker_instance : int
            unique to each worker, to ensure independent seed
        """
        if not self.seeded:
            self.randseed = self.parameters["GlobalSeed"] * (worker_instance+1)
            self.rng = np.random.default_rng(self.randseed)
            self.rng_int = self.rng.integers(low=100, high=10000)
            self.seeded=True
    
    def randint(self)->int:
        """Generate random integer.
            Gives exactly the same result each time unless reseed=True

        Returns
        -------
        int
            a random integer
        """
        if not self.seeded:
            print("NOT SEEDED!!")
            exit(-1)
        
        if self.parameters["FreshSeed"]:
            self.rng_int = self.rng.integers(low=100, high=10000)
        return str(self.rng_int)
    
    def expansion(self,T:float)->np.ndarray:
        """Return the relative x,y,z thermal expansion,
            using the data provided in the XML file

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        np.ndarray, shape (3,)
            relative x,y,z thermal expansion
        """
        scale = np.ones(3) 
        scale += self.parameters["LinearThermalExpansion"]*T
        scale += self.parameters["QuadraticThermalExpansion"]*T*T
        return scale
    
    def info(self)->str:
        """Return all parameters as formatted string

        Returns
        -------
        str
            the parameter string
        """
        result = """
            PAFI parameter information:
            
            Axes:"""
        for k,v in self.axes.items():
            result+=f"""
                {k} : 
                    {v}"""
        result += """
        
            Scripts:"""
        for k,v in self.scripts.items():
            result += f"""
                {k} : {self.parse_script(v)}"""
        result += f"""
            Wildcard Potential Type:
                {self.PotentialType}
            
            Wildcard Potential:
                {self.PotentialLocation}
            
            Wildcard Elements:
                {self.Species}

            Pathway:
                Directory:
                    {self.PathwayDirectory}
                Files:"""
        for p in self.PathwayConfigurations:
            result += f"""
                    {p}"""
        
        result += """


            Parameters"""
        for k,v in self.parameters.items():
            result += f"""
                {k} : {v}"""
        
        result += f"""
                seeded : {self.seeded}"""
        result += f"""
                minValidResults : {self.minValidResults}"""
        result += f"""

        """
        return result
    
    def to_dict(self) -> dict:
        """Export axes and parameters as a nested dictionary

        Returns
        -------
        dict
            Dictionary-of-dictionaries with two keys 'axes' and 'parameters'
        """
        return super().to_dict()
    
    

        

            


