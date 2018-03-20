import numpy as np

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

    @property
    def surfaceCutoff(self):
        return self._surfaceCutoff

    @surfaceCutoff.getter
    def surfaceCutoff(self):
        if self._surfaceCutoff is None:
            # insert here clever way to automatically do that
            pass

        return self._surfaceCutoff

    @surfaceCutoff.setter
    def surfaceCutoff(self, value):
        self._surfaceCutoff = value

    def coarseGrainedDensity(self, points):
        """
        calculate the coarse grained density at a set of specified input points
        """

        # this one needs to be highly efficient (borrow old routine from C code)
        # only pass those atoms actually used in surface construction
        # --> create special array for that!? cached_property?

        return None

    def distanceToSurface(self, points):
        """
        calculate the distance to the ILI using equality constraints
        """

        return None

default_zeta = {
    'O' : 2.5, # Angstrom
    'H' : 0.,
    }
