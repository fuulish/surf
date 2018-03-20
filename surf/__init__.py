import numpy as np
from _surf import coarse_grained_density
from ase.units import Bohr

class ILI(object):
    """
    Instant Liquid Interfaces (ILI) Calculator

    Parameters
    ----------
    atoms : ASE Atoms object
        atoms constituting the system
    mask : list of bool
        (sub)set of atoms used in building the ILI
    zeta : list of floats
        width parameters for exponent used in constructing the ILI
    surfaceCutoff : float
        value of half the bulk density
    """

    def __init__(self, atoms, mask=None, zeta=None, surfaceCutoff=None):
        self.atoms = atoms
        self.mask = mask
        self.zeta = zeta
        self._surfaceCutoff = surfaceCutoff

        if self.mask is None:
            self.mask = [True]*len(atoms)

        # NOTE: equivalent integer mask:
        self.imask = [i for i, v in enumerate(self.mask) if v]

        if self.zeta is None:
            self.zeta = np.array([default_zeta[z] for z in self.atoms.get_chemical_symbols()])

        self.zeta = np.array(self.zeta, dtype='float64')

    @property
    def surfaceCutoff(self):
        return self._surfaceCutoff

    @surfaceCutoff.getter
    def surfaceCutoff(self):
        if self._surfaceCutoff is None:

            cgd = self.coarseGrainedDensity(self.atoms[self.imask].positions)
            self._surfaceCutoff = np.average(cgd) / 2.

        return self._surfaceCutoff

    @surfaceCutoff.setter
    def surfaceCutoff(self, value):
        self._surfaceCutoff = value

    def coarseGrainedDensity(self, points):
        """
        calculate the coarse grained density at a set of specified input points
        """

        points = points #.flatten().astype('float64')
        pos = self.atoms[self.imask].positions.flatten().astype('float64')
        zeta = self.zeta[self.imask]
        natoms = len(zeta)
        pbc = np.diag(self.atoms.get_cell()).astype('float64')

        cgd = []

        for point in points:
            mepos = point.flatten().astype('float64')
            grad = np.zeros(3, dtype='float64')

            cgd.append(coarse_grained_density(mepos, pos, zeta, natoms, pbc, 1, grad))

        return np.array(cgd)

    def distanceToSurface(self, points):
        """
        calculate the distance to the ILI using equality constraints
        """

        return None

    @staticmethod
    def pointsFromCube(cubefile):
        xc, yc, zc = cubefile.coord_grid(0.,0.,0.)

        xc = xc[:,0,0].astype('float64')
        yc = yc[0,:,0].astype('float64')
        zc = zc[0,0,:].astype('float64')

        X, Y, Z = np.meshgrid(xc, yc, zc, indexing='ij')

        points = []

        for x, y, z in zip(X.flatten(), Y.flatten(), Z.flatten()):
            points.append(np.array([x,y,z])*Bohr)

        points = np.array(points)

        return points

default_zeta = {
    'O' : 2.5, # Angstrom
    'H' : 0.,
    }
