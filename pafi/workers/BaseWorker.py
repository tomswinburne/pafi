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

        self.kB = 8.617e-5 # where should I put this...

        self.worker_instance = worker_instance
        self.comm = comm
        self.local_rank = comm.Get_rank()
        self.roots = roots
        self.rank = rank
        self.parameters = parameters
        self.error_count = 0
        self.scale = np.eye(3)
        self.out_width=16
        self.natoms=0
        self.nlocal=0
        self.offset=0
        self.name="BaseWorker"
        self.has_errors = False
        
        self.alpha = 1.0 # mixing of cell and position distance
        

        self.Cell = np.eye(3)
        self.X = np.zeros((1,3))
        self.Volume = 1.0
        self.Periodicity = np.zeros(3,bool)
        self.invCell = np.zeros((3,3))
        self.dCell = np.zeros((3,3))
        self.depsilon = np.zeros((3,3))
        
        self.parameters.seed(worker_instance)
    
    def thermal_expansion(self,T:float=0) -> None:
        """Rescale the LAMMPS supercell *not coordinates*
            according to the provided thermal expansion data

        Parameters
        ----------
        T : float, optional
            temperature in K, by default 0
        """
        s = self.parameters.expansion(T)
        if isinstance(s,float):
            self.scale = np.eye(3) * scale
        elif isinstance(s,np.ndarray):
            assert s.size==3, "expansion must be 3-array or float"
            self.scale = np.diag(s)
        return self.scale

    def load_and_update(self,file_path:os.PathLike[str]) -> np.ndarray:
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
    
    def update(self,C=None,dC=None,X=None,T=None)->None:
        """
            Set state information
        """
        if C is not None:
            self.set_cell(C)
        
        if X is not None:
            self.set_positions(X)
        
        self.Cell, self.Periodicity = self.get_cell()
        self.Volume = np.linalg.det(self.Cell)
        self.invCell = np.linalg.inv(self.Cell)
        
        self.X = self.get_positions()

        self.natoms = self.X.shape[0]
        
        # cell strain (if given)
        if dC is not None:
            self.dCell = dC.copy()
            self.depsilon = self.dCell@self.invCell
            self.norm_dC = np.sqrt((dC**2).sum())
        
        # tangent (if given)
        if T is not None:
            self.norm_T = 1.0 * np.linalg.norm(T)


    def set_path_and_update(self,r:float=0):
        # positions: x and d_u[x,y,z]
        path_S  = self.Spline_S(r,nu=0).reshape((-1,3))
        path_C  = self.Spline_C(r,nu=0).reshape((3,3)) @ self.scale
        path_X  = path_S@path_C
        path_dC = self.Spline_C(r,nu=1).reshape((3,3)) @ self.scale
        
        path_T = self.Spline_S(r,nu=1).reshape((-1,3))@path_C
        path_T -= path_T.mean(0)
        # sets norm_T and norm_dC as well
        self.update(C=path_C,dC=path_dC,X=path_X,T=path_T)
        
        # tangent normalized
        path_T /= self.norm_T

        # dtangent rescaled
        path_dT = self.Spline_S(r,nu=2).reshape((-1,3))@path_C
        path_dT /= self.norm_T**2 + self.norm_dC**2 * self.alpha
        return path_X, path_T, path_dT

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
        self.load_and_update(pc[0])
        invC = self.invCell.copy()
        S = self.X.reshape((-1,3))@self.invCell
        S -= np.floor(S)@np.diag(self.Periodicity)

        all_S = [S.flatten()]
        all_C = [self.Cell.flatten()]
        r_dist = [0.0]
        for p in pc[1:]:
            self.load_and_update(p)
            S = self.X.reshape((-1,3))@self.invCell
            dS = S-all_S[0].reshape((-1,3))
            dS -= np.floor(dS+0.5)@np.diag(self.Periodicity)
            dC = self.Cell@invC - np.eye(3) # engineering distortion
            r_dist += [ np.sqrt( (dS**2).sum() + self.alpha * (dC**2).sum() ) ]
            all_S += [dS.flatten()+all_S[0]]
            all_C += [self.Cell.flatten()]

        # determine distance TODO: symmetric??
        if self.parameters("RealMEPDist"):
            self.r_dist = np.array(r_dist) / r_dist[-1]
        else:
            self.r_dist = np.linspace(0.,1.,all_X.shape[0])

        # splining - thank you scipy for 'axis' !!
        all_S = np.array([S.flatten() for S in all_S])
        all_C = np.array([C.flatten() for C in all_C])

        bc = self.parameters("CubicSplineBoundaryConditions")
        assert bc in ['clamped','not-a-knot','natural']
        self.Spline_S = CubicSpline(self.r_dist,all_S,axis=0,bc_type=bc)
        self.Spline_C = CubicSpline(self.r_dist,all_C,axis=0,bc_type=bc)
        del all_S, all_C # save a bit of memory

    def close(self)->None:
        pass
            


        

