################################################################################
# Carlos Leon
# Wilfrid Laurier University
################################################################################

import os
import numpy as np

################################################################################
# taken from https://wiki.fysik.dtu.dk/ase/_modules/ase/io/vasp.html#write_vasp
def _symbol_count_from_symbols(symbols):
    """Reduce list of chemical symbols into compact VASP notation

    args:
        symbols (iterable of str)

    returns:
        list of pairs [(el1, c1), (el2, c2), ...]
    """
    sc = []
    psym = symbols[0]
    count = 0
    for sym in symbols:
        if sym != psym:
            sc.append((psym, count))
            psym = sym
            count = 1
        else:
            count += 1
    sc.append((psym, count))
    return sc

def getVec(matrix):
	lista = []
	for vec in matrix:
		for e in vec:
			lista = lista + [e]
	#
	# lista = list(map( str, lista ))
	return lista

def getTypes(symbols):
	sc = _symbol_count_from_symbols(symbols)
	i = 0
	types = []
	for _, count in sc:
		for j in range(count):
			types = types + [i]
		#
		i += 1
	#
	types = list(map( str, types ))
	return types
#

# improved from https://wiki.fysik.dtu.dk/ase/_modules/ase/build/tools.html#sort
def getSortedIndices(atoms, tags=None):
    """Return a new Atoms object with sorted atomic order. The default
    is to order according to chemical symbols, but if *tags* is not
    None, it will be used instead. A stable sorting algorithm is used.

    Example:

    >>> from ase.build import bulk
    >>> # Two unit cells of NaCl:
    >>> a = 5.64
    >>> nacl = bulk('NaCl', 'rocksalt', a=a) * (2, 1, 1)
    >>> nacl.get_chemical_symbols()
    ['Na', 'Cl', 'Na', 'Cl']
    >>> nacl_sorted = sort(nacl)
    >>> nacl_sorted.get_chemical_symbols()
    ['Cl', 'Cl', 'Na', 'Na']
    >>> np.all(nacl_sorted.cell == nacl.cell)
    True
    """
    if tags is None:
        tags = atoms.get_chemical_symbols()
    else:
        tags = list(tags)
    #
    mydict = {'C':1000, 'O':1001, 'B':1002, 'N':1003}
    tags = [ mydict[t] for t in tags ]
    #    
    deco = sorted([(tag, i) for i, tag in enumerate(tags)])
    indices = [i for tag, i in deco]
    # return atoms[indices]
    return indices
#

################################################################################
#%%

def atoms2cfg(atoms):
    """
    Save atoms ASE object into a cfg MTP format file.
    """
    cell  = atoms.get_cell() # I had to do this to avoid QE approximations due to conversion from cart to alat units and back
    pos   = atoms.get_positions()

    #!##############################################################################
    #! ******* VERY IMPORTANT!!! ********
    #! sort alphabetically, so I would obtain C first, H as second type, and O as third
    #! it seems you will always want to put 3 ATOMS types (C,H,O) to get appropriate diff.cfg !!! for training
    #! put 1 Oxygen far away from the ribbon in case you don't want to consider it ;) and check the periodicity so that O is always far away from all nanoribbons periodically repeated
    #! from ase.build import sort
    #! atoms = sort(atoms)
    #! I won't use the previous 2 lines because the atoms.get_forces() are lost once they get sorted. this is because the calculator?? is abandoned when doing the slicing, I suppose.
    #! Solution: use get the indices that sort the atoms:
    indices = getSortedIndices(atoms)
    #!##############################################################################

    nAtoms = len(atoms)
    cell   = getVec( cell )
    pos    = getVec( pos[indices] )
    types  = getTypes( np.array(atoms.get_chemical_symbols())[indices] )
    ##

    t = "   "
    # pos = list(map( str, pos ))
    # forces = list(map( str, forces ))
    s  = "BEGIN_CFG\n Size\n"
    s += "    " + str(nAtoms) + "\n"
    s += " SuperCell\n"

    for i in range(3):
        s += "         " + str(cell[0 + (i * 3)]) + t + str(cell[1 + (i * 3)]) + t + str(cell[2 + (i * 3)]) + "\n"
    #
    s += " AtomData:  id type	  cartes_x	cartes_y      cartes_z\n"

    for i in range(nAtoms):
        s += "             " + str(i + 1) + t + types[i]
        s += '%16.7f %16.7f %16.7f\n' % (pos[0 + (i*3)], pos[1 + (i*3)], pos[2 + (i*3)] )
    #
    s += "END_CFG\n\n"
    #
    return s
#



def atoms2cfg2(atoms, energy):
    """
    Save atoms ASE object into a cfg MTP format file.
    """
    cell  = atoms.get_cell() # I had to do this to avoid QE approximations due to conversion from cart to alat units and back
    pos   = atoms.get_positions()

    #!##############################################################################
    #! ******* VERY IMPORTANT!!! ********
    #! sort alphabetically, so I would obtain C first, H as second type, and O as third
    #! it seems you will always want to put 3 ATOMS types (C,H,O) to get appropriate diff.cfg !!! for training
    #! put 1 Oxygen far away from the ribbon in case you don't want to consider it ;) and check the periodicity so that O is always far away from all nanoribbons periodically repeated
    #! from ase.build import sort
    #! atoms = sort(atoms)
    #! I won't use the previous 2 lines because the atoms.get_forces() are lost once they get sorted. this is because the calculator?? is abandoned when doing the slicing, I suppose.
    #! Solution: use get the indices that sort the atoms:
    indices = getSortedIndices(atoms)
    #!##############################################################################

    nAtoms = len(atoms)
    cell   = getVec( cell )
    pos    = getVec( pos[indices] )
    types  = getTypes( np.array(atoms.get_chemical_symbols())[indices] )
    ##

    t = "   "
    # pos = list(map( str, pos ))
    # forces = list(map( str, forces ))
    s  = "BEGIN_CFG\n Size\n"
    s += "    " + str(nAtoms) + "\n"
    s += " SuperCell\n"

    for i in range(3):
        s += "         " + str(cell[0 + (i * 3)]) + t + str(cell[1 + (i * 3)]) + t + str(cell[2 + (i * 3)]) + "\n"
    #
    s += " AtomData:  id type	  cartes_x	cartes_y      cartes_z     fx fy fz\n"

    for i in range(nAtoms):
        s += "             " + str(i + 1) + t + types[i]
        s += '%16.7f %16.7f %16.7f\n' % (pos[0 + (i*3)], pos[1 + (i*3)], pos[2 + (i*3)] )
        s += '%16.7f %16.7f %16.7f\n' % ( 0.0, 0.0, 0.0 )
    #
    s += " Energy\n"
    s += "        " + str(energy) + "\n"
    s += "END_CFG\n\n"
    #
    return s
#




def atoms2cfg(atoms):
    """
    Save atoms ASE object into a cfg MTP format file.
    """
    cell  = atoms.get_cell() # I had to do this to avoid QE approximations due to conversion from cart to alat units and back
    pos   = atoms.get_positions()

    #!##############################################################################
    #! ******* VERY IMPORTANT!!! ********
    #! sort alphabetically, so I would obtain C first, H as second type, and O as third
    #! it seems you will always want to put 3 ATOMS types (C,H,O) to get appropriate diff.cfg !!! for training
    #! put 1 Oxygen far away from the ribbon in case you don't want to consider it ;) and check the periodicity so that O is always far away from all nanoribbons periodically repeated
    #! from ase.build import sort
    #! atoms = sort(atoms)
    #! I won't use the previous 2 lines because the atoms.get_forces() are lost once they get sorted. this is because the calculator?? is abandoned when doing the slicing, I suppose.
    #! Solution: use get the indices that sort the atoms:
    indices = getSortedIndices(atoms)
    #!##############################################################################

    nAtoms = len(atoms)
    cell   = getVec( cell )
    pos    = getVec( pos[indices] )
    types  = getTypes( np.array(atoms.get_chemical_symbols())[indices] )
    ##

    t = "   "
    # pos = list(map( str, pos ))
    # forces = list(map( str, forces ))
    s  = "BEGIN_CFG\n Size\n"
    s += "    " + str(nAtoms) + "\n"
    s += " SuperCell\n"

    for i in range(3):
        s += "         " + str(cell[0 + (i * 3)]) + t + str(cell[1 + (i * 3)]) + t + str(cell[2 + (i * 3)]) + "\n"
    #
    s += " AtomData:  id type	  cartes_x	cartes_y      cartes_z\n"

    for i in range(nAtoms):
        s += "             " + str(i + 1) + t + types[i]
        s += '%16.7f %16.7f %16.7f\n' % (pos[0 + (i*3)], pos[1 + (i*3)], pos[2 + (i*3)] )
    #
    s += "END_CFG\n\n"
    #
    return s
#



def dict2cfg(d, dictOfSpecies):
    """
    Save atoms dictionary into a cfg MTP format file.
    """
    nAtoms = d["size"]
    c      = d["cell"]
    p      = d["positions"]
    f      = d["outputs"]["forces"]
    energy = d["outputs"]["energy"]
    st = d["outputs"]["virial_stress"]


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

    s += " AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n"

    for i in range(nAtoms):
        s += "             " + str(i + 1) + t + types[i]
        s += '  %16.7f %16.7f %16.7f     ' % (p[i][0], p[i][1], p[i][2] )
        s += '  %16.7f %16.7f %16.7f\n' % (f[i][0], f[i][1], f[i][2] )
    #
    s += " Energy\n"
    s += "        " + str(energy) + "\n"
    s += " PlusStress:  xx          yy          zz          yz          xz          xy\n"
    s += '              %16.7f %16.7f %16.7f %16.7f %16.7f %16.7f\n'    % (st[0], st[1], st[2], st[3], st[4], st[5])
    s += "END_CFG\n\n"
    #
    return s
#



