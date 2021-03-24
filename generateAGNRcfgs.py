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

def getCfgs(atoms, mindist, nCfgs, stdev):
    atoms_original = copy.deepcopy(atoms)
    cfgString = ""
    countCfgs = 0
    for i in range(1000): # try 1000 times, until reaching the desired number of configurations `nCfgs`
        atoms.rattle(stdev, seed=i) # see rattle() in https://wiki.fysik.dtu.dk/ase/ase/atoms.html
        if getMinDist(atoms) > mindist: # disregard those structures with atoms too close
            cfgString += atoms2cfg(atoms) 
            countCfgs += 1
            if countCfgs == nCfgs:
                atoms = copy.deepcopy(atoms_original)
                break
            #
        #
        atoms = copy.deepcopy(atoms_original)
        #
        # the following lines does not get random values for x, y, z 
        # newPositions = perturbAllPositions(atoms.get_positions(), amplitud)
        # atoms.set_positions(newPositions)
        # cfgString += atoms2cfg(atoms)
        # atoms = copy.deepcopy(atoms_original)
    #
    print("I generated " + str(countCfgs) + " cfgs.")
    f = open("to_relax.cfg", "w")
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

def getCfgsWithRandAtom(atoms, atom, mindist, nCfgs, stdev):
    atoms_original = copy.deepcopy(atoms)
    pos_original   = copy.deepcopy(atom.position)
    cfgString = ""
    countCfgs = 0
    for i in range(1000): # try 1000 times, until reaching the desired number of configurations `nCfgs`
        _rattleAtom(atom, stdev, seed=i)
        atoms.append(atom)
        if getMinDist(atoms) > mindist: # disregard those structures with atoms too close
            cfgString += atoms2cfg(atoms) 
            countCfgs += 1
            if countCfgs == nCfgs:
                atoms = copy.deepcopy(atoms_original)
                atom.position = copy.deepcopy(pos_original)
                break
            #
        #
        atoms = copy.deepcopy(atoms_original)
        atom.position = copy.deepcopy(pos_original)
    #
    #     atom.position = perturbPosition(atom.position, amplitud)
    #     atoms.append(atom)
    #     cfgString += atoms2cfg(atoms)
    #     atoms = copy.deepcopy(atoms_original)
    #     atom.position = pos_original
    # #
    print("I generated " + str(countCfgs) + " cfgs.")
    f = open("to_relax.cfg", "a")
    f.write(cfgString)
    f.close()
#

def getMinDist(atoms):
    distMatrix = atoms.get_all_distances()
    distances = [e for e in distMatrix.flatten() if e != 0.0]
    mindist1 = np.amin(distances)
    #
    distMatrix = atoms.get_all_distances(mic=True)
    distances = [e for e in distMatrix.flatten() if e != 0.0]
    mindist2 = np.amin(distances)
    #
    return min(mindist1, mindist2)
#

#%%
from ase.build import graphene_nanoribbon

mindist = 0.5
vacuum = 9.0
stdev  = 0.5
nCfgs  = 50
#
atoms = graphene_nanoribbon(2, 1, type='armchair', saturated=False, C_H=1.1, C_C=1.4, vacuum=vacuum,  magnetic=True, initial_mag=1.12)
getCfgs(atoms, mindist=mindist, nCfgs=nCfgs, stdev=stdev)

#%%
# add oxygen at random positions over the unit cell space
from ase import Atom

r = atoms.get_positions()
rhalf = (r[1] + r[2]) / 2

atom = Atom("O", position=rhalf )
getCfgsWithRandAtom(atoms, atom, mindist=mindist, nCfgs=nCfgs, stdev=vacuum)


#%%
# Visualization
from ase.build import graphene_nanoribbon
atoms = graphene_nanoribbon(2, 3, type='armchair', saturated=True, C_H=1.1, C_C=1.4, vacuum=vacuum,  magnetic=True, initial_mag=1.12)
import nglview as nv
# nv.show_ase(atoms)
v = nv.show_ase(atoms)
v.background = 'black'
v


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
dictOfSpecies = {"0":"C", "1":"O"}
dict = read_cfgs("relaxed.cfg", dictOfSpecies)
# %%
from ase import Atoms
d = dict[0]
species = d['species']
cell = d['cell']
positions = d['positions']

atoms = Atoms("".join(species), cell=cell, positions=positions, pbc=(1,0,0))
# %%
import nglview as nv
v = nv.show_ase(atoms)
v.background = 'black'
v
# %%
from ase.io.espresso import write_espresso_in
QEin = "temporal.in"
f = open(QEin, "w")
write_espresso_in(fd=f, atoms=atoms, input_data=None, pseudopotentials=None, kspacing=None, kpts=None, koffset=(0, 0, 0), crystal_coordinates=False)
f.close()

# %%
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
fig, ax = plt.subplots()
plot_atoms(atoms, ax, radii=0.3, rotation=('90x,0y,0z'))
fig.savefig("ase_atoms.png")
#%%

# %%
from ase import Atoms
d = dict[0]
species = d['species']
cell = d['cell']
positions = d['positions']
atoms = Atoms("".join(species), cell=cell, positions=positions, pbc=(1,0,0))
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
fig, ax = plt.subplots()
plot_atoms(atoms, ax, radii=0.3, rotation=('90x,0y,0z'))
# fig.savefig("ase_atoms.png")

# %%
