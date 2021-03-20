#%%
import copy
from random import random
import numpy as np
# seed random number generator
from aseAtoms2cfg import atoms2cfg

def getAmplifiedRand(amplitud=1):
    return amplitud * ( random() - 0.5 )
#
def perturbPosition(position, amplitud):
    position += getAmplifiedRand(amplitud)
    return position
#
def perturbAllPositions(positions, amplitud):
    for r in positions:
        r = perturbPosition(r, amplitud)
    #
    return positions
#

def getCfgs(atoms, nCfgs, stdev):
    atoms_original = copy.deepcopy(atoms)
    cfgString = ""
    for i in range(nCfgs):
        atoms.rattle(stdev, seed=i) # see rattle() in https://wiki.fysik.dtu.dk/ase/ase/atoms.html
        cfgString += atoms2cfg(atoms) 
        atoms = copy.deepcopy(atoms_original)
        #
        # the following lines does not get random values for x, y, z 
        # newPositions = perturbAllPositions(atoms.get_positions(), amplitud)
        # atoms.set_positions(newPositions)
        # cfgString += atoms2cfg(atoms)
        # atoms = copy.deepcopy(atoms_original)
    #
    f = open("toRelax.cfg", "w")
    f.write(cfgString)
    f.close()
#
def _rattleAtom(atom, stdev=0.001, seed=None, rng=None):
    """Randomly displace position of an atom of type `Atom`.

    This method adds a random displacement to the atomic position,
    taking a possible constraint into account??.  The random numbers are
    drawn from a normal distribution of standard deviation stdev.

    For a parallel calculation, it is important to use the same
    seed on all processors!  """

    if seed is not None and rng is not None:
        raise ValueError('Please do not provide both seed and rng.')

    if rng is None:
        if seed is None:
            seed = 42
        rng = np.random.RandomState(seed)
    pos = atom.position
    atom.position = pos + rng.normal(scale=stdev, size=3)
#

def getCfgsWithRandAtom(atoms, atom, nCfgs, stdev):
    atoms_original = copy.deepcopy(atoms)
    pos_original   = copy.deepcopy(atom.position)
    cfgString = ""
    for i in range(nCfgs):
        _rattleAtom(atom, stdev, seed=i)
        atoms.append(atom)
        cfgString += atoms2cfg(atoms) 
        atoms = copy.deepcopy(atoms_original)
        atom.position = copy.deepcopy(pos_original)


    #     atom.position = perturbPosition(atom.position, amplitud)
    #     atoms.append(atom)
    #     cfgString += atoms2cfg(atoms)
    #     atoms = copy.deepcopy(atoms_original)
    #     atom.position = pos_original
    # #
    f = open("toRelax.cfg", "a")
    f.write(cfgString)
    f.close()
#   
#%%
from ase.build import graphene_nanoribbon

vacuum = 6.0
atoms = graphene_nanoribbon(1, 1, type='armchair', saturated=False, C_H=1.1, C_C=1.4, vacuum=vacuum,  magnetic=True, initial_mag=1.12)
#
getCfgs(atoms, nCfgs=50, stdev=1.0)

#%%
# add oxygen at random positions over the unit cell space
from ase import Atom

r = atoms.get_positions()
rhalf = (r[1] + r[2]) / 2

atom = Atom("O", position=rhalf )
getCfgsWithRandAtom(atoms, atom, nCfgs=50, stdev=vacuum)



# print(atoms.get_positions())
# print("")
# newPositions = perturbAllPositions(atoms.get_positions(), 1.0)
# print(newPositions)
# print("")
# atoms.set_positions(newPositions)
# print(atoms.get_positions())




# %%
import cfg2ase
import importlib
importlib.reload(cfg2ase)

from cfg2ase import read_cfgs
dict = read_cfgs("toRelax.cfg")
# %%
