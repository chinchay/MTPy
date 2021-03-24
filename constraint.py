import numpy as np
class MyConstraint:
    """Constrain an atom to move along a given direction only."""
    def __init__(self, a, direction):
        self.a = a
        self.dir = direction / sqrt(np.dot(direction, direction))

    def adjust_positions(self, atoms, newpositions):
        step = newpositions[self.a] - atoms.positions[self.a]
        step = np.dot(step, self.dir)
        newpositions[self.a] = atoms.positions[self.a] + step * self.dir

    def adjust_forces(self, atoms, forces):
        forces[self.a] = self.dir * np.dot(forces[self.a], self.dir)