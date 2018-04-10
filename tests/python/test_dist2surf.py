import numpy as np
from ase.io import read
from ase.units import Bohr

from surf import ILI


def test_dist2surf(shared_datadir):

    atoms = read(shared_datadir / 'conf.xyz')
    atoms.set_pbc([True, True, True])
    atoms.set_cell([36, 36, 100])

    mask = [atom.symbol == 'O' for atom in atoms]
    ili = ILI(atoms, mask=mask, density_cutoff=0.002223/Bohr**3)

    dsts = ili.distanceToPoint(atoms[mask].positions)
    dsts_gsl = ili.distanceToPoint(atoms[mask].positions, gsl=True)

    np.testing.assert_allclose(dsts, dsts_gsl, rtol=20)

    dsts_ref = np.loadtxt(str(shared_datadir / 'refdist.dat'))
    dsts_ref = dsts_ref[:, 1] * Bohr

    np.testing.assert_allclose(dsts, dsts_ref, rtol=5)
    np.testing.assert_allclose(dsts_gsl, dsts_ref, rtol=5)
