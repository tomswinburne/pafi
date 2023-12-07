import numpy as np
import os
from typing import List
from mpi4py import MPI
from ..parsers.PAFIParser import PAFIParser
from scipy.interpolate import CubicSpline

class BaseWorker:
    """Basic PAFI Worker

        Parameters
        ----------
        comm : MPI.Intracomm
            MPI communicator
        parameters : PAFIParser
            Predefined or custom  PAFIParser object
        worker_instance : int
            unique worker rank
        rank : int
            global rank (for MPI safety)
        roots : List[int]
            list of master ranks  (for MPI safety)
        """
    def __init__(self, comm : MPI.Intracomm,
                 parameters:PAFIParser,
                 worker_instance:int,
                 rank:int,
                 roots:List[int]) -> None:
        self.worker_instance = worker_instance
        self.comm = comm
        self.local_rank = comm.Get_rank()
        self.roots = roots
        self.rank = rank
        self.parameters = parameters
        self.error_count = 0
        self.scale = np.ones(3)
        self.out_width=16
        self.natoms=0
        self.nlocal=0
        self.offset=0
        self.x=None
        self.name="BaseWorker"
        self.has_errors = False
        self.has_cell_data=False
        self.Cell = None
        self.Periodicity = None
        self.invCell = None
        self.parameters.seed(worker_instance)
    
    def load_config(self,file_path:os.PathLike[str]) -> np.ndarray:
        """Placeholder function to load in file and return configuration
        Overwritten by LAMMPSWorker
        
        Parameters
        ----------
        file_path : os.PathLike[str]
            path to file

        Returns
        -------
        np.ndarray, shape (natoms,3)
            configuration vector
        """
        assert os.path.exists(file_path)
        return np.loadtxt(file_path)

    def pbc(self,X:np.ndarray,central:bool=True)->np.ndarray:
        """Minimum image convention, using cell data
            
        Parameters
        ----------
        X : np.ndarray, shape (natoms,3)
            configuration vector
        central : bool, optional
            map scaled coordinates to [-.5,.5] if True, else [0,1], by default True

        Returns
        -------
        np.ndarray
            wrapped vector
        """
        if not self.has_cell_data:
            return X
        else:
            sX = X.reshape((-1,3))@self.invCell
            sX -= np.floor(sX+0.5*int(central))@np.diag(self.Periodicity)
            return (sX@self.Cell).reshape((X.shape))
        
        return X
    def pbc_dist(self,X:np.ndarray,axis:None|int=None)->float|np.ndarray:
        """Minimum image distance
        
        Parameters
        ----------
        X : np.ndarray, shape (natoms,3)
            configuration vector
        axis : None | int, optional
            as in np.linalg.norm, by default None
        
        Returns
        -------
        float|np.ndarray
            norm of the vector or vector(s)
        """
        return np.linalg.norm(self.pbc(X),axis=axis)
    
    def make_path(self):
        """Make the splined PAFI path via scipy.interpolate.CubicSpline

            All parameters are read in from XML file

            Lots to do here in principle, TODO: parallel i/o and splining
            Only really a problem with memory limitations, say 2GB / core.
            This implies 300M double coordinates == 1M atoms, 100 planes.
            Typical large-scale use - 150k atoms, 20 planes.
            So memory-heavy but nothing problematic so far, leaving for future
        """

        # load configurations
        pc = self.parameters.PathwayConfigurations
        all_X = [self.pbc(self.load_config(pc[0]),central=False)]
        for p in pc[1:]:
            all_X += [self.pbc(self.load_config(p)-all_X[0])+all_X[0]]
            
        # determine distance TODO: symmetric??
        if self.parameters("RealMEPDist"):
            self.r_dist = np.array([self.pbc_dist(X-all_X[0]) for X in all_X])
            self.r_dist /= self.r_dist[-1]
        else:
            self.r_dist = np.linspace(0.,1.,all_X.shape[0])

        # splining - thank you scipy for 'axis' !!
        all_X = np.array([X.flatten() for X in all_X]) # shape=(nknots,natoms)
        bc = self.parameters("CubicSplineBoundaryConditions")
        assert bc in ['clamped','not-a-knot','natural']

        self.Spline_X = CubicSpline(self.r_dist,all_X,axis=0,bc_type=bc)
        del all_X # save a bit of memory

    def pathway(self,r:float,nu:int=0,
                scale:float|np.ndarray[3]=1.0)->np.ndarray:
        """Evaluate the PAFI pathway

        Parameters
        ----------
        r : float
            reaction coordinate, should be in [0,1]
        nu : int, optional
            derivative order, by default 0
        scale : float | np.ndarray[3], optional
            thermal expansion, by default 1.0

        Returns
        -------
        np.ndarray, shape (natoms,3)
            pathway configuration
        """
        if isinstance(scale,float):
            scale = np.eye(3)*float(scale)
        else:
            scale = np.diag(scale)
        return self.Spline_X(r,nu=nu).reshape((-1,3))@scale
    
    def close(self)->None:
        pass
            


        

