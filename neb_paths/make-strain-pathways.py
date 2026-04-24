# %%
from ase import Atoms
from ase.io import write,read 
from ase.lattice.cubic import BodyCenteredCubicFactory
from ase.filters import StrainFilter
from ase.calculators.eam import EAM
from ase.calculators.lammpslib import LAMMPSlib
import numpy as np
import os
import matplotlib.pyplot as plt 
import nglview as nv


a_bcc = 3.179040316319168  # Typical BCC Fe lattice constant in Å
a_fcc =  3.84249488769937


def EquilibriumCubic(Factory,weights,a0=3.6,size=[3,3,3],element='W'):
        LML_pot = LMLPotential(element=element, coefficients=weights)
        calc = LML_pot.make_lammpslib_calc()
        atoms = Cartesian(a0,[1,1,1],element,Factory)
        atoms.calc = calc 
        
        stress_ratio = 1.0/np.abs(calc.get_stress()[:3].mean())

        
        sf = StrainFilter(atoms)
        opt = PreconLBFGS(sf, precon=None, logfile=None)
        opt.run(fmax=1e-4, smax=1e-3)  # max force in eV/A

        stress_ratio *= np.abs(calc.get_stress()[:3].mean())
        
        lattice_constant = atoms.cell.lengths()[0]
        
        atoms = Cartesian(lattice_constant,size,element,Factory)
        atoms.calc = calc
        N = atoms.get_positions().shape[0]        
        print("VOL",np.linalg.det(atoms.cell)/N,N)
        
        return lattice_constant, np.round(stress_ratio,3)

def generate_bain_path(a_bcc=3.18, a_fcc = 3.85, n_cell = 2, n_images=7, traj='w-bain-path.traj'):
    F = BodyCenteredCubicFactory()
    atoms_bcc = F(directions=[[1,0,0], [0,1,0], [0,0,1]], 
                    size=[n_cell]*3,
                    symbol='W', pbc=(1,1,1),
                    latticeconstant=a_bcc)

    c_over_a = np.linspace(1.0, np.sqrt(2), n_images)
    
    vol_per_atom = np.linspace(1.0, (a_fcc**3 / 4.0)/(a_bcc**3 / 2.0), n_images)

    length_per_atom = vol_per_atom ** (1.0/3.0)

    V0 = atoms_bcc.get_volume()

    atoms_list = []

    for i, c in enumerate(c_over_a):
    
        # Create a new atoms object for each image
        atoms = atoms_bcc.copy()
        # Set the c/a ratio
        cell = atoms.get_cell()
        V0 = np.linalg.det(cell)
        cell[2, 2] = c* cell[0, 0]
        cell *= (V0 / np.linalg.det(cell)) ** (1.0/3.0)
        cell *= length_per_atom[i]
        atoms.set_cell(cell, scale_atoms=True)
        atoms_list.append(atoms)
    
    write(traj, atoms_list)
    

generate_bain_path(a_bcc=a_bcc, a_fcc=a_fcc, n_cell=2, n_images=17)


def generate_strain_path(a_bcc=3.18, a_fcc = 3.85, n_cell = 2, n_images=7, traj='w-strain-path.traj'):
    F = BodyCenteredCubicFactory()
    atoms_bcc = F(directions=[[1,0,0], [0,1,0], [0,0,1]], 
                    size=[n_cell]*3,
                    symbol='W', pbc=(1,1,1),
                    latticeconstant=a_bcc)

    e_array = np.linspace(-0.0,0.07,n_images)
    atoms_list = []
    for i, e in enumerate(e_array):
    
        # Create a new atoms object for each image
        atoms = atoms_bcc.copy()

    
        # Set the c/a ratio
        cell = np.array(atoms.get_cell())
        cell[0][1] = e*cell[0][0]
        
        atoms.set_cell(cell, scale_atoms=True)

        atoms_list.append(atoms)

    write(traj, atoms_list)

generate_shear_path(a_bcc=a_bcc, a_fcc=a_fcc, n_cell=2, n_images=17)


def generate_expand_path(a_bcc=3.18, a_fcc = 3.85, n_cell = 2, n_images=7, traj='w-expand-path.traj'):
    F = BodyCenteredCubicFactory()
    atoms_bcc = F(directions=[[1,0,0], [0,1,0], [0,0,1]], 
                    size=[n_cell]*3,
                    symbol='W', pbc=(1,1,1),
                    latticeconstant=a_bcc)

    e_array = 1.0+np.linspace(-0.03,0.07,n_images)


    atoms_list = []

    for i, e in enumerate(e_array):
    
        # Create a new atoms object for each image
        atoms = atoms_bcc.copy()

    
        # Set the c/a ratio
        cell = atoms.get_cell()
        cell *= e
        atoms.set_cell(cell, scale_atoms=True)

        atoms_list.append(atoms)
    
        
    write(traj, atoms_list)

generate_expand_path(a_bcc=a_bcc, a_fcc=a_fcc, n_cell=2, n_images=17)


