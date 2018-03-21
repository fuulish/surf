#!/usr/bin/env python

import numpy as np
from surf import ILI
from ase.io import read, write
from ase.units import Bohr
from ase import Atoms, Atom

from cube import Cube

atoms = read('conf.xyz')
atoms.set_pbc([True, True, True])
atoms.set_cell([36, 36, 100])

mask = [atom.symbol == 'O' for atom in atoms]
surface = ILI(atoms, mask=mask, surfaceCutoff=0.002223/Bohr**3)

dsts = surface.distanceToSurface(atoms[mask].positions)

np.savetxt('distances.txt', dsts)

prev = np.loadtxt('refdist.dat')
prev = prev[:,1] * Bohr

np.testing.assert_allclose(dsts, prev, rtol=5)
