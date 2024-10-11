import unittest,os
SELF_DIR = os.path.dirname(os.path.abspath(__file__))
DUMP_DIR = os.path.join(SELF_DIR,"./")

class PAFITestCase(unittest.TestCase):
	def test_import(self):
		# Add a name to the test case
		"""
		Test case for importing necessary modules
		"""
		import numpy # just to check import
		import scipy # just to check import
		import pandas # just to check import
		from mpi4py import MPI
		from pafi.workers.LAMMPSWorker import wrappedlammps
		self.assertTrue(MPI)
		self.assertTrue(wrappedlammps)
		comm = MPI.COMM_WORLD
		rank = comm.Get_rank()
		size = comm.Get_size()
		
		lmp = wrappedlammps(comm=comm,\
					 cmdargs=["-log","none","-screen","none"])
		lmp.close()
		
		self.assertTrue(rank >= 0)
		self.assertTrue(size > 0)
		print(""" 

		""")

	def test_load_parameters(self,nowrite=True):
		# Add a name to the test case
		"""
		Test case for loading parameters
		"""
		from pafi.parsers.PAFIParser import PAFIParser
		self.assertTrue(PAFIParser)
		self.parameters = PAFIParser(postprocessing=nowrite)
		self.parameters.set("DumpFolder",DUMP_DIR)
		# set path wildcard for Fe SIA
		self.parameters.set_pathway(\
			os.path.join(SELF_DIR,"EAM-SIA-Fe/image_*.dat"))
		# set EAM interatomic potential
		self.parameters.set_potential(\
			os.path.join(SELF_DIR,"EAM-SIA-Fe/Fe.eam.fs"))
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
		print(""" 




		""")
	
	def test_zero_temperature_run(self):
		# Add a name to the test case
		"""
		Test case for running at zero temperature
		"""
		import numpy as np
		from mpi4py import MPI
		from pafi.managers.PAFIManager import PAFIManager
		self.assertTrue(PAFIManager)
		self.test_load_parameters(False)
		manager = PAFIManager(world=MPI.COMM_WORLD,parameters=self.parameters)
		results = manager.run()
		manager.close()

		test_gradient = np.array([0.,0.614145,0.846192,0.59226,0.089912,
							 -0.426347,-0.853405, -0.792713, 0.])


		sim_gradient = np.array(results['FreeEnergyGradient'])
		
		self.assertTrue(np.isclose(sim_gradient,test_gradient).all())
		print(""" 




		""")

if __name__ == '__main__':
	unittest.main()
