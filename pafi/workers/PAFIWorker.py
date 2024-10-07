import numpy as np
import os
from mpi4py import MPI
from typing import Any, List
from ..parsers.PAFIParser import PAFIParser
from .LAMMPSWorker import LAMMPSWorker
from ..results.ResultsHolder import ResultsHolder

class PAFIWorker(LAMMPSWorker):
    """
        The hyperplane-constrained sampling run.
        results: ResultsHolder instance
            essentially a dictionary, defining both 
            input and output data- temeperature, reaction coordinate
            and other are defined, and others are filled.
            *Crucially, if a parameter is specified in results, 
            it overrides any value in the parameter file*
            Must have at least `ReactionCoordinate` and `Temperature` defined
        
        The PAFI workflow can be summarised as:
            1) Execute `PreRun` script
            2) Apply `fix_pafi` constraint at defined `ReactionCoordinate`
            3) Execute `PreTherm` script
            4) Thermalization for `ThermSteps` steps at `Temperature`
            5) Execute `constrained_average()` function
                In standard PAFI this runs for `SampleSteps` steps and 
                time averages the output of `fix_pafi`, as shown below.
                See https://docs.lammps.org/fix_pafi.html for details.
            6) Extract the average displacment from path if `PostDump==1`
            7) Minimize in-plane to check system returns to path.
                The check is the max per-atom displacement : `MaxJump`
                If `MaxJump` is larger than `MaxJumpMaxJumpThresh` then 
                the sample is retained but marked as `Valid=False`
            8) Execute `PostRun` script
        
        The `contrained_average()` function is therefore suitable for
        any form of hyperplane-constrained averaging.
        """
    
    def __init__(self, comm: MPI.Intracomm, 
                 parameters: PAFIParser, tag: int,
                 rank: int, roots: List[int]) -> None:
        super().__init__(comm, parameters, tag, rank, roots)
        
    
    def constrained_average(self,results:ResultsHolder)->ResultsHolder:
        """
        
        Parameters
        ----------
        results: ResultsHolder instance
            used to contain results and custom input paramaters
            we want these inputs to override self.parameters()
        
        Returns
        ----------
        ResultsHolder instance
            Returns the input data and all output data appended 
            as dictionary key,value pairs
        """
        # ensure results inputs should override self.parameters()
        parameters = lambda k: results(k)\
            if results.has_key(k) else self.parameters(k)
        fixname = self.setup_pafi_average(parameters("SampleSteps"),"avepafi")
        self.run_commands("run %d" % parameters("SampleSteps"))
        results = self.extract_pafi_data(results,fixname)
        
        return results
    
    def sample(self,results:ResultsHolder)->ResultsHolder:
        """
        Main sampling run.
        
        1) Execute `PreRun` script
        2) Apply `fix_pafi` constraint at defined `ReactionCoordinate`
        3) Execute `PreTherm` script
        4) Thermalization for `ThermSteps` steps at `Temperature`
        
        5) Execute some averaging function
            In standard PAFI this runs for `SampleSteps` steps and 
            time averages the output of `fix_pafi`, as shown below.
            See https://docs.lammps.org/fix_pafi.html for details.
        
        6) Extract the average displacment from path if `PostDump==1`
        7) Minimize in-plane to check system returns to path.
            The check is the max per-atom displacement : `MaxJump`
            If `MaxJump` is larger than `MaxJumpMaxJumpThresh` then 
            the sample is retained but marked as `Valid=False`
        8) Execute `PostRun` script

        
        Parameters
        ----------
        results: ResultsHolder instance
            used to contain results and custom input paramaters
            we want these inputs to override self.parameters()
        
        Returns
        ----------
        ResultsHolder instance
            Returns the input data and all output data appended 
            as dictionary key,value pairs
        """
        
    
        results = self.standard_pafi_pre_average(results)
        results = self.constrained_average(results)
        results = self.standard_pafi_post_average(results)

        return results
    

    def setup_pafi_average(self,ave_steps:int,fixname="avepafi")->str:
        """Helper function to establish PAFI average
        Parameters
        ----------
        
        ave_steps : int
            steps for ave/time
        fixname:str, optional
            the fix name, default "avepafi"
        Returns
        ---------
        str:
            the fix name
        """
        self.run_commands(f"""
            fix {fixname} all ave/time 1 {ave_steps} {ave_steps} f_pafi[*]
        """)
        return fixname
    
    def extract_pafi_data(self,results:ResultsHolder,
                          name:str="avepafi")->ResultsHolder:
        """Helper functions to extract PAFI data

        Parameters
        ----------
        results : ResultsHolder instance
            add data and returns
        name : str, optional
            fix name, default "avepafi"

        Returns
        -------
        ResultsHolder instance
        """
        res = {}
        fix_data = self.extract_fix(name,size=4)
        res['FreeEnergyGradient'] = -fix_data[0] * self.norm_t
        res['FreeEnergyGradientVariance'] = fix_data[1]**2 * self.norm_t**2 - res['FreeEnergyGradient']**2
        res['avePsi'] = fix_data[2]
        res['dXTangent'] = fix_data[3]
        results.set_dict(res)
        self.run_commands(f"unfix {name}")
        return results
    
    
    def standard_pafi_pre_average(self,results:ResultsHolder)->ResultsHolder:
        """Helper functions for standard PAFI

        Parameters
        ----------
        results : ResultsHolder instance
            custom input data overrides parameters
        
        Returns
        ----------
        results : ResultsHolder instance
            add data and returns custom inputs as well    
        """
        assert results.has_key("ReactionCoordinate")
        assert results.has_key("Temperature")
        
        self.kB = 8.617e-5 # where should I put this...

        # we want any input paramaters in results() to override parameters()
        parameters = lambda k: results(k) \
            if results.has_key(k) else self.parameters(k)
        
        # initialize on hyperplane
        r = results("ReactionCoordinate")
        T = results("Temperature")
        
        self.initialize_hyperplane(r,0.)
        # TODO: check order in public PAFI
        self.run_script("PreRun",results)
        self.initialize_hyperplane(r,T)
        n_atoms = self.get_natoms()

        # the PAFI fix
        gamma = parameters("Friction")
        overdamped = parameters("OverDamped")
        seed = self.parameters.randint()
        cmd = f"fix pafi all pafi __pafipath {T} {gamma} {seed} "
        cmd += f"overdamped {overdamped} com 1\nrun 0"
        self.run_commands(cmd)

        # pre minimize (optional)
        if parameters("PreMin"):
            min_steps = parameters("MinSteps")
            self.run_commands(f"""
                min_style fire
                minimize 0 0.0001 {min_steps} {min_steps}
            """)
        if parameters("PostMin"):
            self.change_x = -self.gather("x",1,3)
        
        # PreThermalize (optional)
        self.run_script("PreTherm",results)
        results.set("MinEnergy",self.get_energy())
        
        # establish temperature time average and thermalize
        pre_therm_hp = self.extract_fix("pafi",size=4)
        steps = parameters("ThermSteps")
        ave_steps = parameters("ThermWindow")
        f_T = "c_thermo_pe" if overdamped==1 else "c_thermo_temp"
        self.run_commands(f"""
            reset_timestep 0
            fix __ae all ave/time 1 {ave_steps} {steps} {f_T}
            run {steps}
        """)
        sampleT = self.extract_fix("__ae")
        
        if overdamped==1:
            sampleT = (sampleT-results("MinEnergy"))/1.5/n_atoms/self.kB
        results.set("preTemperature",sampleT)
        self.run_commands(f"""
            unfix __ae
            run 0""")
        
        # main sampling run
        steps = parameters("SampleSteps")
        self.run_commands(f"""
            reset_timestep 0
            fix __ae all ave/time 1 {steps} {steps} {f_T}
            """)
        if parameters("PostDump"):
            self.run_commands(f"""
            fix pafiax all ave/atom 1 {steps} {steps} x y z
            """)
        return results
    
    def standard_pafi_post_average(self,results:ResultsHolder)->ResultsHolder:
        """Helper functions for standard PAFI

        Parameters
        ----------
        results : ResultsHolder instance
            custom input data overrides parameters
        
        Returns
        ----------
        results : ResultsHolder instance
            add data and returns custom inputs as well    
        """
        # we want any input paramaters in results() to override self.parameters()
        parameters = lambda k: results(k) if results.has_key(k) else self.parameters(k)
        n_atoms = self.get_natoms()
        
        # get final temperature
        sampleT = self.extract_fix("__ae")
        if parameters("OverDamped")==1:
            sampleT = (sampleT-results("MinEnergy"))/1.5/n_atoms/self.kB
        results.set("postTemperature",sampleT)
        self.run_commands(f"""
            unfix __ae
            run 0""")
        
        # average positions
        if parameters("PostDump"):
            dx = self.gather("f_pafiax",type=1,count=3)
            dx[:,0] -= self.gather("f_ux",type=1,count=1).flatten()
            dx[:,1] -= self.gather("f_uy",type=1,count=1).flatten()
            dx[:,2] -= self.gather("f_uz",type=1,count=1).flatten()
            dx = self.pbc(dx)
            results.set("MaxDev",np.linalg.norm(dx,axis=1).max())
            results.set("Dev",dx.copy())
            del dx
            self.run_commands("unfix pafiax")
        
        # minimize to test for MaxJump
        if parameters("PostMin"):
            min_steps = parameters("MinSteps")
        else:
            min_steps = 1
        self.run_commands(f"""
            min_style fire
            minimize 0 0.0001 {min_steps} {min_steps}
        """)
        self.change_x += self.gather("x",1,3)
        results.set("MaxJump",self.pbc_dist(self.change_x,axis=1).max())
        results.set("Valid",bool(results("MaxJump")<parameters("MaxJumpThresh")))
        del self.change_x
        # unfix hyperplane
        self.run_commands("unfix pafi")
        self.run_script("PostRun",results)
        # rescale back.... not sure if this is required
        r = results("ReactionCoordinate")
        self.initialize_hyperplane(r,0.0)
        return results
  
            



        





        


