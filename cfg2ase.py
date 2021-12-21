#***********************************************************************
#* Program:
#*    cfg2ase.py
#* Author:
#*    Carlos Leon, Wilfrid Laurier University
#* Summary:
#*    This is a library to help with the reading of cfg format file (from
#*    MTP) to ASE atoms object.
#*************************************************************************/

from ase import Atom
import numpy as np
import re

def formatify(string):
    return [float(s) for s in string.split()]
#

def getSize(block):
    size_pattern = re.compile("Size\n(.*?)\n SuperCell", re.S | re.I)
    size_str = size_pattern.findall(block)[0]
    return int(size_str.lstrip())
#

def getCell(block):
    cell_pattern = re.compile("SuperCell\n(.*?)\n AtomData", re.S | re.I)
    cell_str = cell_pattern.findall(block)[0]
    return list(map(formatify, cell_str.split("\n")))
#

def getSpeciesPositionsForces(block, dictOfSpecies):
    """
    species, positions, forces = getSpeciesPositionsForces(block)
    """
    position_pattern = re.compile("fz\n(.*?)\n Energy", re.S)
    matrix_str = position_pattern.findall(block)[0]

    matrix_strLines = matrix_str.split("\n")
    species = [ dictOfSpecies[ line.split()[1] ] for line in matrix_strLines ]

    matrix_floats = np.array(list(map(formatify, matrix_str.split("\n"))))    
    positions = matrix_floats[:, 2:5]
    forces    = matrix_floats[:, 5:8].tolist()
    
    # if len(position_pattern.findall(block)) > 0:
        # matrix_str = position_pattern.findall(block)[0]
        # matrix_floats = np.array(list(map(formatify, matrix_str.split("\n"))))
        # species   = np.array(self.elements)[matrix_floats[:, 1].astype(np.int64)]
        # positions = matrix_floats[:, 2:5]
        # forces    = matrix_floats[:, 5:8].tolist()
    # else:
    #     position_pattern = re.compile("cartes_z\n(.*?)\nEND_CFG", re.S)
    #     print(position_pattern.findall(block))
    #     matrix_str = position_pattern.findall(block)[0]
        
    #     matrix_floats = np.array(list(map(formatify, line_str.split("\n"))))
    #     species   = np.array(self.elements)[matrix_floats[:, 1].astype(np.int64)]
    #     positions = matrix_floats[:, 2:5]
    #     forces    = []
    # #
    return species, positions, forces
#

def getEnergy(block):
    energy_pattern = re.compile("Energy\n(.*?)\n (?=PlusStress|Stress)", re.S)
    energy_str = energy_pattern.findall(block)[0]
    return float(energy_str.lstrip())
#

def getStress(block):
    stress_pattern = re.compile("xy\n(.*?)(?=\n|$)", re.S)
    stress_str = stress_pattern.findall(block)[0]
    virial_stress = formatify(stress_str)
    return virial_stress
#

# inspired on ttps://github.com/materialsvirtuallab/maml/blob/master/maml/apps/pes/_mtp.py
def read_cfgs(filename, dictOfSpecies):
    """
    Args:
        filename (str): The configuration file to be read.
    """
    with open(filename, "r") as f:
        # lines = f.readlines() # << it won't work in the `for` loop :/
        letters = f.read() # 
    #
    block_pattern = re.compile("BEGIN_CFG\n(.*?)\nEND_CFG", re.S)
    
    data_pool = []
    for block in block_pattern.findall(letters):
        d = {"outputs": {}}
        d["size"] = getSize(block)
        d["cell"] = getCell(block)
        species, positions, forces = getSpeciesPositionsForces(block, dictOfSpecies)
        d["species"]   = species
        d["positions"] = positions
        assert d["size"] == len(species)
        #
        d["outputs"]["energy"] = getEnergy(block)
        d["outputs"]["forces"] = forces
        d["outputs"]["virial_stress"] = getStress(block)
        #
        data_pool.append(d)
    #
    return data_pool
#

def species2elements(species, dictionary):
    elements = copy.deepcopy(species)
    for (i, s) in enumerate(species):
        elements[i] = dictionary[s]
    #
    return elements
#

def getSpeciesPositions(block, dictOfSpecies):
    """
    species, positions = getSpeciesPositions(block)
    """

    position_pattern1 = re.compile("cartes_z\n(.*?)\n Feature", re.S)
    position_pattern2 = re.compile("cartes_z\n(.*?)(?=$)", re.S) #  `$` means end of string

    vec1 = position_pattern1.findall(block)
    vec2 = position_pattern2.findall(block)

    vec = vec1 if len(vec1) > 0 else vec2

    matrix_str = vec[0]

    matrix_strLines = matrix_str.split("\n")
    species = [ dictOfSpecies[ line.split()[1] ] for line in matrix_strLines ]

    matrix_floats = np.array(list(map(formatify, matrix_str.split("\n"))))    
    positions = matrix_floats[:, 2:5]

    return species, positions
#


# inspired on ttps://github.com/materialsvirtuallab/maml/blob/master/maml/apps/pes/_mtp.py
def read_cfgs_asinput(filename, dictOfSpecies):
    """
    Args:
        filename (str): The configuration file to be read.
    """
    with open(filename, "r") as f:
        # lines = f.readlines() # << it won't work in the `for` loop :/
        letters = f.read() # 
    #
    block_pattern = re.compile("BEGIN_CFG\n(.*?)\nEND_CFG", re.S)
    
    data_pool = []
    for block in block_pattern.findall(letters):
        d = {}
        d["size"] = getSize(block)
        d["cell"] = getCell(block)
        species, positions = getSpeciesPositions(block, dictOfSpecies)
        d["species"]   = species
        d["positions"] = positions
        assert d["size"] == len(species)
        #
        data_pool.append(d)
    #
    return data_pool
#

def read_cfgs_general(filename, dictOfSpecies):
    """
    Args:
        filename (str): The configuration file to be read.
    """
    with open(filename, "r") as f:
        # lines = f.readlines() # << it won't work in the `for` loop :/
        letters = f.read() # 
    #
    block_pattern = re.compile("BEGIN_CFG\n(.*?)\nEND_CFG", re.S)
    
    data_pool = []
    for block in block_pattern.findall(letters):
        d = {}
        d["size"] = getSize(block)
        d["cell"] = getCell(block)

        try:
            species, positions, forces = getSpeciesPositionsForces(block, dictOfSpecies)
            d["species"]   = species
            d["positions"] = positions
            assert d["size"] == len(species)
            #
            d["outputs"] = {}
            d["outputs"]["energy"] = getEnergy(block)
            d["outputs"]["forces"] = forces
            d["outputs"]["virial_stress"] = getStress(block)
        except:
            species, positions = getSpeciesPositions(block, dictOfSpecies)
            d["species"]   = species
            d["positions"] = positions
            assert d["size"] == len(species)
            # 
        # 
        data_pool.append(d)
    #
    return data_pool
#

#%%
from ase import Atoms, Atom

def dict2atoms(dictCfg, pbc):
    species = dictCfg['species']
    cell = dictCfg['cell']
    positions = dictCfg['positions']
    atoms = Atoms("".join(species), cell=cell, positions=positions, pbc=pbc)
    return atoms
#


def dict2cfg(d, dictOfSpecies):
    """
    Save atoms dictionary into a cfg MTP format file.
    """
    nAtoms = d["size"]
    c      = d["cell"]
    p      = d["positions"]
    # f      = d["outputs"]["forces"]
    # energy = d["outputs"]["energy"]
    # st = d["outputs"]["virial_stress"]


    dictOfTypes = {}
    for key in dictOfSpecies:
        val = dictOfSpecies[key]
        dictOfTypes[val] = key
    #

    types = []
    for specie in d["species"]:
        types.append(dictOfTypes[specie])
    #
    
    ################################################

    t  = "   "
    t2 = "         "
    # pos = list(map( str, pos ))
    # forces = list(map( str, forces ))
    s  = "BEGIN_CFG\n Size\n"
    s += "    " + str(nAtoms) + "\n"
    s += " SuperCell\n"

    s += '         %16.7f %16.7f %16.7f\n' % (c[0][0], c[0][1], c[0][2] )
    s += '         %16.7f %16.7f %16.7f\n' % (c[1][0], c[1][1], c[1][2] )
    s += '         %16.7f %16.7f %16.7f\n' % (c[2][0], c[2][1], c[2][2] )

    # s += " AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n"
    s += " AtomData:  id type       cartes_x      cartes_y      cartes_z\n"

    for i in range(nAtoms):
        s += "             " + str(i + 1) + t + types[i]
        s += '  %16.7f %16.7f %16.7f     ' % (p[i][0], p[i][1], p[i][2] )
        s += "\n"
        # s += '  %16.7f %16.7f %16.7f\n' % (f[i][0], f[i][1], f[i][2] )
    #
    # s += " Energy\n"
    # s += "        " + str(energy) + "\n"
    # s += " PlusStress:  xx          yy          zz          yz          xz          xy\n"
    # s += '              %16.7f %16.7f %16.7f %16.7f %16.7f %16.7f\n'    % (st[0], st[1], st[2], st[3], st[4], st[5])
    s += "END_CFG\n\n"
    #
    return s
#




#%%
def getDx(atoms, i, j):
    return (atoms[i].position - atoms[j].position)[0]
#

def belongs(x, x1, x2):
    return ( (x1 < x) and (x < x2) )
#

def sacudeParalelo(atoms, i, j, eps, a):
    dx = getDx(atoms, i, j)
    if not belongs(dx, -eps, eps):
        atoms[j].position[0] += a
        dx = getDx(atoms, i, j)
        if not belongs(dx, -eps, eps):
            atoms[j].position[0] -= (2*a)
#

def zigzagear(atoms):
    a, b, alat, _, _, _ = atoms.cell.cellpar()
    eps = 0.2

    sacudeParalelo(atoms, 1, 2, eps, a)
    sacudeParalelo(atoms, 1, 5, eps, a)
    sacudeParalelo(atoms, 1, 6, eps, a)
    sacudeParalelo(atoms, 1, 9, eps, a)
    sacudeParalelo(atoms, 1, 10, eps, a)
    sacudeParalelo(atoms, 1, 13, eps, a)
    sacudeParalelo(atoms, 1, 14, eps, a)

    sacudeParalelo(atoms, 0, 3, eps, a)
    sacudeParalelo(atoms, 0, 4, eps, a)
    sacudeParalelo(atoms, 0, 7, eps, a)
    sacudeParalelo(atoms, 0, 8, eps, a)
    sacudeParalelo(atoms, 0, 11, eps, a)
    sacudeParalelo(atoms, 0, 12, eps, a)
    sacudeParalelo(atoms, 0, 15, eps, a)

    # Now oxygens:
    sacudeParalelo(atoms, 0, 16, eps, a)
    sacudeParalelo(atoms, 0, 17, eps, a)
#

def mywrap_C802(cfgFileName, dictOfSpecies, fileNameOut):
    dict = read_cfgs_general(cfgFileName, dictOfSpecies)
    cfgString = ""
    for d in dict:
        atoms = dict2atoms(d)

        # Wrap positions to unit cell.
        # See https://wiki.fysik.dtu.dk/ase/_modules/ase/atoms.html#Atoms.set_scaled_positions
        # See ase.geometry.wrap_positions() in https://wiki.fysik.dtu.dk/ase/ase/geometry.html
        atoms.wrap(pbc=[True, True, True]) # <<-- MTP and QE will receive more work or even will give erroneous mindist (mlp mindist) or bad convergence (QE)

        zigzagear(atoms)

        d["positions"] = atoms.positions
        cfgString += dict2cfg(d, dictOfSpecies)
    #
    f = open(fileNameOut, "w")
    f.write(cfgString)
    f.close()    
#    



#