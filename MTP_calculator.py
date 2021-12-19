# inspired on https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/lj.html#LennardJones
import numpy as np

from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
# from ase.stress import full_3x3_to_voigt_6_stress


class MTP_calculator(Calculator):
    """
    Calculator for MTP potential
    """
    implemented_properties = ['energy']
    # implemented_properties = ['energy', 'energies', 'forces', 'free_energy']
    # implemented_properties += ['stress', 'stresses']  # bulk properties

    def __init__(self, filename):
        """
        Parameters
        ----------
        filename: string
            File with coefficients/parameters
        """
        Calculator.__init__(self, filename)
    #

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        # if properties is None:
        #     properties = self.implemented_properties
        # #
        # Calculator.calculate(self, atoms, properties, system_changes)

        self.results['energy'] = 101.0
    #
#
