import numpy as np
import os
from typing import Any, List
from ..parsers.PAFIParser import PAFIParser
from mpi4py import MPI
from lammps import lammps,LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,LMP_TYPE_SCALAR
from .BaseWorker import BaseWorker
from ..results.ResultsHolder import ResultsHolder

class wrappedlammps(lammps):
    """
        Wrapper for lammps for safe MPI initialization

        Try to pass MPI communicator to lammps, if fails 
        we assume lammps is compiled without MPI support


    """
    def __init__(self,name='',cmdargs=None,ptr=None,comm=None):
        try:
            super().__init__(name=name,cmdargs=cmdargs,comm=comm)
            self.has_mpi_comm = True
        except Exception as ae:
            # in case we have conda-lammps installed 
            super().__init__(name=name,cmdargs=cmdargs,comm=None)
            self.has_mpi_comm = False
            # check if we are only assigning one core if no MPI support
            assert not self.lib.lammps_config_has_mpi_support()
            if comm is not None:
                assert comm.Get_size()==1



class LAMMPSWorker(BaseWorker):
    """LAMMPS worker for PAFI, inheriting BaseWorker

        Initialization routine:
        - runs "Input" script, loading first configuration
        - extracts cell data
        - defines various functions wrapping computes etc.
        - defines initialize_hyperplane function to set plane

        Parameters
        ----------
        comm : MPI.Intracomm
            MPI communicator
        parameters : PAFIParser
            Predefined or custom PAFIParser object
        worker_instance : int
            unique worker rank
        """
    def __init__(self, comm: MPI.Intracomm, 
                 parameters: PAFIParser, worker_instance: int,
                 rank: int, roots: List[int]) -> None:
        super().__init__(comm, parameters, worker_instance, rank, roots)
        
        self.name = "LAMMPSWorker"
        self.last_error_message = ""
        self.start_lammps()

        if (not self.L.has_mpi_comm) and comm.Get_size()>1:
            self.has_errors = True
            print("LAMMPS HAS NO MPI SUPPORT- ONE CORE/WORKER ONLY!")
            return
        
        if self.has_errors:
            print("ERROR STARTING LAMMPS!",self.last_error_message)
            return
        start_config = self.parameters.PathwayConfigurations[0]
        
        # TODO abstract
        self.run_script("Input")
        if self.has_errors:
            print("ERROR RUNNING INPUT SCRIPT!",self.last_error_message)
            return
        self.update()
        
        # See initialize_hyperplane
        self.scale = np.ones(3)
        self.made_fix=False
        self.made_compute=False
        self.make_path()
    
        
    def start_lammps(self)->None:
        """Initialize LAMMPS instance

            Optionally 

        """
        if self.parameters("LogLammps"):
            logfile = 'log.lammps.%d' % self.worker_instance 
        else:
            logfile = 'none'
        try:
            cmdargs = ['-screen','none','-log',logfile]
            self.L = wrappedlammps(comm=self.comm,cmdargs=cmdargs)
            self.check_lammps_compatibility()
        except Exception as ae:
            print("Couldn't load LAMMPS!",ae)
            self.has_errors = True
    
    def check_lammps_compatibility(self)->None:
        """
            Ensure LAMMPS is new enough and has fix_pafi
        """
        lammps_release_int = self.L.version()
        if lammps_release_int<20201101: 
            if self.local_rank == 0:
                print("Require LAMMPS version > 01Nov2020!")
            return
        pafi_package = "EXTRA-FIX" if lammps_release_int>=20210728 else "USER-MISC"
        self.has_pafi = self.L.has_package(pafi_package)
        if not self.has_pafi and self.local_rank==0:
            print("Cannot find PAFI package in LAMMPS!")
            self.has_errors = True
        
    def run_script(self,key:str,arguments:None|dict|ResultsHolder=None)->None:
        """Run a script defined in the XML
            Important to replace any %wildcards% if they are there!
            
        Parameters
        ----------
        key : str
            script key from XML
        arguments : None | dict | ResultsHolder, optional
            will be used to replace wildcards, by default None
        """
        if key in self.parameters.scripts:
            if self.parameters("Verbose")>0 and self.rank==0:
                print(f"RUNNING SCRIPT {key}")
        
            script = self.parameters.parse_script(key,arguments=arguments)
            self.run_commands(script)

    def run_commands(self,cmds : str | List[str]) -> bool:
        """
            Run LAMMPS commands line by line, checking for errors
        """
        cmd_list = cmds.splitlines() if isinstance(cmds,str) else cmds
        for cmd in cmd_list:
            try:
                if self.parameters("Verbose")>0 and self.rank==0:
                    print(f"TRYING COMMAND {cmd}")
                self.L.command(cmd)
            except Exception as ae:
                if self.local_rank==0:
                    message = f"LAMMPS ERROR: {cmd} {ae}"
                else:
                    message = None
                self.last_error_message = ae
                raise SyntaxError(message)
    
    def gather(self,name:str,type:None|int=None,count:None|int=None)->np.ndarray:
        """Wrapper of LAMMPS gather()
           

        Parameters
        ----------
        name : str
            name of data
        type : None | int, optional
            type of array, 0:integer or 1:double, by default None.
            If None, an attempt will be made to determine autonomously.
        count : None | int, optional
            number of data per atom, by default None. 
            If None, an attempt will be made to determine autonomously.

        Returns
        -------
        np.ndarray
            the LAMMPS data

        Raises
        ------
        ValueError
            if name not found
        """
        name=name.lower()
        if name in ['x','f','v'] or 'f_' in name:
            if type is None:
                type = 1
            if count is None:
                count = 3
        elif name in ['id','type','image']:
            if type is None:
                type = 0
            if count is None:
                count = 1
        if type is None or count is None:
            raise ValueError("Error in gather: type or count is None")
        
        try:
            res = self.L.gather(name,type,count)
        except Exception as ae:
            if self.local_rank==0:
                print("Error in gather:",ae)
            self.last_error_message = ae
        return np.ctypeslib.as_array(res).reshape((-1,count))
    
    def scatter(self,name:str,data:np.ndarray)->None:
        """Scatter data to LAMMPS
            Assume ordered with ID

        Parameters
        ----------
        name : str
            name of array
        data : np.ndarray
            numpy array of data. Will be flattened.
        """
        if np.issubdtype(data.dtype,int):
            type = 0
        elif np.issubdtype(data.dtype,float):
            type = 1
        count = data.shape[1] if len(data.shape)>1 else 1
        try:
            if self.parameters("Verbose")>0 and self.rank==0:
                print("Scattering",name,type,count,data.shape)
            self.L.scatter(name,type,count,
                           np.ctypeslib.as_ctypes(data.flatten()))
        except Exception as ae:
            if self.local_rank==0:
                print("Error in scatter:",ae)
            self.last_error_message = ae
    
    def get_positions(self):
        """
        Get positions
        """
        return self.gather("x",1,3)
    
    def set_positions(self,x):
        """
        Set positions
        """
        self.scatter("x",x)
        return None

    def get_cell(self):
        """
        Get Cell
        """
        boxlo,boxhi,xy,yz,xz,pbc,box_change = self.L.extract_box()
        Periodicity = np.array([bool(pbc[i]) for i in range(3)],bool)
        Cell = np.zeros((3,3))
        for cell_j in range(3):
            Cell[cell_j][cell_j] = boxhi[cell_j]-boxlo[cell_j]
        Cell[0][1] = xy
        Cell[0][2] = xz
        Cell[1][2] = yz
        return Cell,Periodicity

    def set_cell(self,C):
        """ 
        Set Cell
        """
        self.run_commands(f"""
                change_box all triclinic

                change_box all x final 0.0 {C[0][0]} y final 0.0 {C[1][1]} z final 0.0 {C[2][2]} xy final {C[0][1]} xz final {C[0][2]} yz final {C[1][2]}

                run 0
            """)
        return None

    def get_natoms(self)->int:
        """Get the atom count

        Returns
        -------
        int
            the atom count
        """
        return self.L.get_natoms()

    def setup_stress_average(self, ave_steps:int)->None:
        """
        
        Helper function to setup stress average
        
        Parameters
        ----------
        ave_steps : int
            steps for ave/time
        
        ASSUMES METAL UNITS!
            Stress is in bar = 1e5 J/m^3
                             = 1e5 (1e19/1.6)eV / (1e30 A^3)
                             = (1.0/1.6) * 1e-6 eV / A^3
                             = ConvFactor eV/A^3
            Vol is in A^3 
            
            => Stress * Vol * ConvFactor is in eV
        
        Returns
        ---------
        None
        
        """
        # assumes metal units !!
        conv = (1/1.6) * (1e-6) 
        self.run_commands(f"""
            variable conv equal {conv}
            variable pxx equal pxx*vol*v_conv
            variable pyy equal pyy*vol*v_conv
            variable pzz equal pzz*vol*v_conv
            variable pxy equal pxy*vol*v_conv
            variable pxz equal pxz*vol*v_conv
            variable pyz equal pyz*vol*v_conv
            fix axx all ave/time 1 {ave_steps} {ave_steps} v_pxx
            fix ayy all ave/time 1 {ave_steps} {ave_steps} v_pyy
            fix azz all ave/time 1 {ave_steps} {ave_steps} v_pzz
            fix axy all ave/time 1 {ave_steps} {ave_steps} v_pxy
            fix axz all ave/time 1 {ave_steps} {ave_steps} v_pxz
            fix ayz all ave/time 1 {ave_steps} {ave_steps} v_pyz
        """)

    def extract_stress_average(self)->np.ndarray:
        """
            See setup_stress_average
        """
        stress_fixes = ['axx', 'ayy', 'azz', 'axy', 'axz', 'ayz']
        stress = np.zeros(6)
        for i, var in enumerate(stress_fixes):
            stress[i] = -self.extract_fix(var,size=1)
        return stress
    
    def unset_stress_average(self)->None:
        self.run_commands(f"""
            unfix axx
            unfix ayy
            unfix azz
            unfix axy
            unfix axz
            unfix ayz
        """)    
    
    def load_and_update(self,file_path:os.PathLike[str])->None:
        """Load a LAMMPS data file with read_data() and return a numpy array of positions

        Parameters
        ----------
        file_path : os.PathLike[str]
            Path to the LAMMPS .dat file.

        Returns
        -------
        np.ndarray, shape (N,3)
            the positions
        """
        self.made_fix = False
        self.make_compute = False
        self.run_commands(f"""
            delete_atoms group all
            read_data {file_path} add merge
        """)
        self.update()

    def extract_compute(self,id:str,vector:bool=True)->float|np.ndarray:
        """    Extract compute from LAMMPS
            
            id : str
        Parameters
        ----------
        id : str
            compute id
        vector : bool, optional
            is the return value a vector, by default True

        Returns
        -------
        np.ndarray or float
           return data
        """
        style = LMP_STYLE_GLOBAL
        type = LMP_TYPE_VECTOR if vector else LMP_TYPE_SCALAR
        assert hasattr(self.L,"numpy")
        try:
            res = self.L.numpy.extract_compute(id,style,type) 
            return np.array(res)
        except Exception as e:
            if self.local_rank==0:
                print("FAIL EXTRACT COMPUTE",e)
            self.close()
            return None
    
    def extract_fix(self,id:str,size:int=1)->float|np.ndarray:
        """Extract fix from LAMMPS and return a numpy array

        Parameters
        ----------
        id : str
            name of fix
        size : int, optional
            number of global entries, by default 1
        Returns
        -------
        float or np.ndarray
            numpy array of data of shape (size,), or float if size=1
        """
        
        style = LMP_STYLE_GLOBAL
        type = LMP_TYPE_VECTOR if size>1 else LMP_TYPE_SCALAR
        assert hasattr(self.L,"numpy")
        try:
            res = lambda i: self.L.numpy.extract_fix(id,style,type,ncol=i)
            if size>1:
                return np.array([res(i) for i in range(size)])
            else:
                return res(0)
        except Exception as e:
            if self.local_rank==0:
                print("FAIL EXTRACT FIX",id,e)
            self.close()
            return None
    
    def get_energy(self)->float:
        """
        Extract the potential energy
        
        Returns
        -------
        float
            The potential energy (pe)
        """
        return self.extract_compute("thermo_pe",vector=False)

    def initialize_hyperplane(self,r:float,T:float)->None:
        """Establish worker on a given plane with the 
            reference pathway via `fix_pafi`. 
            Also computes the tangent magnitude.
        Parameters
        ----------
        r : float
            The reaction coordinate, typically confined to [0,1]
        T : float
            The temperature in K. This will apply thermal expansion to the pathway.
        """
        
        # check for __pafipath fix to store the path data
        # (path : d_u[x,y,z],tangent: d_n[x,y,z],dtangent: d_dn[x,y,z])
        
        if not self.made_fix:
            self.run_commands(f"""
                fix __pafipath all property/atom d_ux d_uy d_uz d_nx d_ny d_nz d_dnx d_dny d_dnz
                run 0
            """)
            self.made_fix=True
        
        self.scale = self.thermal_expansion(T) # updates self.scale 
        
        path_X,path_T,path_dT = self.set_path_and_update(r)
        
        # fill path
        for i,c in enumerate(["d_ux","d_uy","d_uz"]):
            self.scatter(c,path_X[:,i])
        
        for i,c in enumerate(["d_nx","d_ny","d_nz"]):
            self.scatter(c,path_T[:,i])
        
        for i,c in enumerate(["d_dnx","d_dny","d_dnz"]):
            self.scatter(c,path_dT[:,i])
        
        #del path_X,path_T,path_dT
        
        self.run_commands("run 0")


        # check for __pafipath compute to make this data accessible
        if not self.made_compute:
            self.run_commands(f"""
            compute __pafipath all property/atom d_ux d_uy d_uz d_nx d_ny d_nz d_dnx d_dny d_dnz
            run 0
            """)
            self.made_compute=True
    
    def close(self) -> None:
        """
            Close down. TODO Memory management??
        """
        super().close()
        self.L.close()

    