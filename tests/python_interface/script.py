#!/usr/bin/env python

import numpy as np
from surf import ILI
from ase.io import read

atoms = read('conf.xyz')
atoms.set_pbc([True, True, True])
atoms.set_cell([36, 36, 100])

mask = [atom.symbol == 'O' for atom in atoms]
surface = ILI(atoms, mask=mask)

points = np.zeros((3,1))
points = atoms.positions.copy()

cgd = surface.coarseGrainedDensity(points)

print 'coarse grained density is: ', cgd
