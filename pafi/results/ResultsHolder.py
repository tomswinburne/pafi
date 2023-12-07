from typing import Any,List

class ResultsHolder:
    """PAFI ResultsHolder object, 
    a dictionary which contains input and output values

    Methods
    -------
    __call__
    has_key
    """
    def __init__(self) -> None:
        self.data = {}
    def items(self)-> Any: # TODO what is the type here?
        """return items for iteration
        Returns
        -------
        dictionary.items()
        """
        return self.data.items()
    def __call__(self,key:str)->Any:
        """_summary_

        Parameters
        ----------
        key : str
            dictionary key

        Returns
        -------
        Any
           dictionary value
        
        Raises
        ------
        ValueError
            If key not found
        """
        if self.has_key(key):
            return self.data[key]
        else:
            raise ValueError(f"No key {key} in Results!")
        
    def set(self,key:str,value:Any)->None:
        """Set key,value pair in internal data

        Parameters
        ----------
        key : str
        value : Any
        """
        self.data[key] = value

    def has_key(self,key:str)->bool:
        """Check if key exists

        Parameters
        ----------
        key : str
           key query

        Returns
        -------
        bool
            True if exists
        """
        return bool(key in self.data.keys())
    
    def set_dict(self,dict:dict)->None:
        """Write or overwrite to results from dictionary
        writes key,value pairs via set()

        Parameters
        ----------
        dict : dict
        """
        for key,value in dict.items():
            self.set(key,value)
    def get_dict(self,keys:List[str],blanks:None|str=None)->dict:
        """Returns dictionary with given keys 

        Parameters
        ----------
        keys : List[str]
            list of keys
        blanks : None or str
            If not None, return str for missing key, default None

        Returns
        -------
        dict
            dictionary with values determined from internal data
        """
        res = {}
        for key in keys:
            if not blanks is None:
                res[key] = self.data[key] if self.has_key(key) else blanks
            else:
                res[key] = self.__call__(key)
        return res