import unittest
import glob
import numpy as np
from pafi import PAFIParser, ResultsProcessor
from mpi4py import MPI

class PAFITestCase(unittest.TestCase):
    def test_load_parameters(self,nowrite=True):
        self.parameters = PAFIParser(postprocessing=nowrite)
        # set path wildcard for Fe SIA
        self.parameters.set_pathway(\
            "../examples/systems/EAM-SIA-Fe/image_*.dat")
        # set EAM interatomic potential
        self.parameters.set_potential(\
            "../examples/systems/EAM-SIA-Fe/Fe.eam.fs")
        # set potential species
        self.parameters.set_species(["Fe"])
        # restrict to zero temperature for testing 
        # TODO: finite temperature tests!
        self.parameters.axes["Temperature"] = [0.]
        # single force call for thermalization and sampling at T=0K
        self.parameters.set("nRepeats",1)
        self.parameters.set("Verbose",0)
        self.parameters.set("SampleSteps",1)
        self.parameters.set("ThermSteps",1)
        self.parameters.set("ThermWindow",1)
        self.assertTrue(self.parameters.ready())
        print(self.parameters.info())
    
    def test_zero_temperature_run(self):
        self.test_load_parameters(False)
        from pafi import PAFIManager
        manager = PAFIManager(world=MPI.COMM_WORLD,parameters=self.parameters)
        results = manager.run()
        manager.close()

        test_gradient = np.array([0.,0.614145,0.846192,0.59226,0.089912,
                             -0.426347,-0.853405, -0.792713, 0.])


        sim_gradient = np.array(results['FreeEnergyGradient'])
        
        self.assertTrue(np.isclose(sim_gradient,test_gradient).all())

if __name__ == '__main__':
    unittest.main()