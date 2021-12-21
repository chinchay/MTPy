import importlib
import utils

# inspired on https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/lj.html#LennardJones
from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
# from ase.stress import full_3x3_to_voigt_6_stress


import numpy as np
from ase import Atoms

importlib.reload(utils)


class MTP_calculator(Calculator):
    """
    Calculator for MTP potential
    """
    implemented_properties = ['energy']
    # implemented_properties = ['energy', 'energies', 'forces', 'free_energy']
    # implemented_properties += ['stress', 'stresses']  # bulk properties

    def __init__(self, filename, dictionaryTypes):
        """
        Parameters
        ----------
        filename: string
            File with coefficients/parameters
        """
        Calculator.__init__(self)
        self.parameters = utils.load(filename)
        self.initializedVecs = utils.init_vecs(self.parameters)
        self.dictionaryTypes = dictionaryTypes
    #

    def calculate(self, atoms, properties=None, system_changes=all_changes):
        # if properties is None:
        #     properties = self.implemented_properties
        # #
        # Calculator.calculate(self, atoms, properties, system_changes)

        # energy = utils.CalcEFS()

        # self.atoms = atoms # the structure `atoms` to wich this calculator is attached
        # print(self.atoms)

        self.atoms = atoms # the structure `atoms` to wich this calculator is attached 
        self.list_max_dist = 0.5 * self.parameters["max_dist"] * np.ones(len(atoms))
        neighborhoods = utils.get_neighborhoods(atoms, self.list_max_dist, self.dictionaryTypes)
        type_centrals = utils.get_type_centrals(atoms, self.dictionaryTypes)
        energy = utils.CalcEFS(self.atoms, neighborhoods, type_centrals, self.parameters, self.initializedVecs)



        
        self.results['energy'] = energy
    #
#
