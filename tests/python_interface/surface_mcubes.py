#!/usr/bin/env python

import mcubes
from ase import Atom
from ase.io import write, read
import numpy as np
from surf import ILI
from ase.units import Bohr

from cube import Cube

atoms = read('conf.xyz')
atoms.set_pbc([True, True, True])
atoms.set_cell([36, 36, 100])

mask = [atom.symbol == 'O' for atom in atoms]
surface = ILI(atoms, mask=mask)

dx = 1.0
lower = (0, 0, 0)
upper = np.diag(atoms.get_cell())
stride = np.array((upper - lower) / dx, dtype='int')

upper = tuple(upper)

f = lambda x, y, z: surface.coarseGrainedDensity([np.array([x,y,z])])
verts, triangles = mcubes.marching_cubes_func(lower, upper, stride[0], stride[1], stride[2], f, surface.surfaceCutoff)

for vert in verts:
    atoms.append(Atom('X', position=vert))
    
write('mcubes_surface.xyz', atoms)
