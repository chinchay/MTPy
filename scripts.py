from random import random
# seed random number generator
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
# def generateNperturbedCfgs(atoms, nCfgs, amplitud):
#     for i in range(nCfgs):
#         newPositions = perturbAllPositions(atoms.get_positions(), amplitud)
#         atoms.set_positions(newPositions)
# #        

# random()generateNperturbedCfgs(atoms, 3, 1.0):

atoms = graphene_nanoribbon(1, 1, type='armchair', saturated=False, C_H=1.1, C_C=1.4, vacuum=6.0,  magnetic=True, initial_mag=1.12)

print(atoms.get_positions())
print("")
newPositions = perturbAllPositions(atoms.get_positions(), 1.0)
print(newPositions)
print("")
atoms.set_positions(newPositions)
print(atoms.get_positions())