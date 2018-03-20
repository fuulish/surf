#!/usr/bin/env python

import numpy as np
from surf import ILI
from ase.io import read
from ase.units import Bohr

from cube import Cube

atoms = read('conf.xyz')
atoms.set_pbc([True, True, True])
atoms.set_cell([36, 36, 100])


mask = [atom.symbol == 'O' for atom in atoms]
surface = ILI(atoms, mask=mask)

points = np.zeros((3,1))
points = atoms.positions.copy()

tcb = Cube('test.cube')
points = surface.pointsFromCube(tcb)

cgd = surface.coarseGrainedDensity(points)

new = tcb.copy()
# new = Cube('empty.cube')

new.data[:] = 0.
new.data[:] = np.reshape(cgd, tcb.data.shape) * Bohr**3

new_data = new.data.flatten()
old_data = tcb.data.flatten()

np.testing.assert_allclose(new_data, old_data, atol=1.e-6)

tcb.write_to_file('output.cube') #, tcb)
