#!/usr/bin/env python

from surf import ILI
from ase.io import read

atoms = read('conf.xyz')
atoms.set_pbc([True, True, True])
atoms.set_cell([36, 36, 100])

mask = [atom.symbol == 'O' for atom in atoms]
surface = ILI(atoms, mask=mask)
cgd = surface.coarseGrainedDensity([[0.,0.,0],])

print 'coarse grained density is: ', cgd
