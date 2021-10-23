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
