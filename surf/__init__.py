import numpy as np
from _surf import coarse_grained_density
from ase.units import Bohr
from scipy.optimize import minimize

try:
    import mcubes
    mcubes_loaded = True
except ImportError:
    mcubes_loaded = False

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
        self.pbc = np.diag(self.atoms.get_cell())

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

    def coarseGrainedDensity(self, points, gradient=False):
        """
        calculate the coarse grained density at a set of specified input points
        """

        points = points #.flatten().astype('float64')
        pos = self.atoms[self.imask].positions.flatten().astype('float64')
        zeta = self.zeta[self.imask]
        natoms = len(zeta)
        pbc = np.diag(self.atoms.get_cell()).astype('float64')

        cgd = np.zeros(len(points))
        grd = np.zeros((len(points), 3))
        calc_grad = [0,1][gradient]

        for i, point in enumerate(points):
            mepos = point.flatten().astype('float64')
            grad = np.zeros(3, dtype='float64')

            # TODO: only calculate gradient if requested
            cgd[i] = coarse_grained_density(mepos, pos, zeta, natoms, pbc, 1, calc_grad, grad)
            grd[i][:] = grad[:]

        if gradient:
            return grd

        return cgd

    def distanceToSurface(self, points):
        """
        calculate the distance to the ILI using equality constraints
        """

        # minimize distance between point and surface
        def func(x, point):
            dst = self.calculate_distance_vector(x, point)
            return (dst**2).sum()

        constraints = ({
            'type' : 'eq',
            'fun' : lambda x: self.coarseGrainedDensity([x]) - self.surfaceCutoff,
        })

        # obtain initial guess from brute-force point search
        sg = self.surfaceGrid(dx=1.)

        dst = np.zeros(len(points))

        for i, point in enumerate(points):

            distToSurf = self.calculate_distance(point, sg)

            x0 = sg[np.argmin(distToSurf)]

            ret = minimize(func, x0, args=(point,), method='SLSQP', constraints=constraints)

            if not ret.success:
                print('Failed to achieve convergence for point #%i: ' %i, point)

            dst[i] = self.calculate_distance(ret.x, point)

        return dst

    def surfaceGrid(self, dx=2., refine=True, marching_cubes=False):
        """
        brute-force estimate of surface position
        """

        if marching_cubes and mcubes_loaded:
            lower = (0, 0, 0)
            upper = self.pbc
            stride = np.array((upper - lower) / dx, dtype='int')
            upper = tuple(upper)

            f = lambda x, y, z: self.coarseGrainedDensity([np.array([x,y,z])])
            verts, triangles = mcubes.marching_cubes_func(lower, upper, \
                               stride[0], stride[1], stride[2], f, self.surfaceCutoff)

            return verts

        # generate grid points

        xg = np.arange(0, self.pbc[0], dx)
        yg = np.arange(0, self.pbc[1], dx)
        zg = np.arange(0, self.pbc[2], dx)

        X, Y = np.meshgrid(xg, yg, indexing='ij')

        surfacePoints = []

        for x, y in zip(X.flatten(), Y.flatten()):
            pts = np.vstack([np.repeat(x, len(zg)), np.repeat(y, len(zg)), zg]).T

            density = self.coarseGrainedDensity(pts)
            gradient = self.coarseGrainedDensity(pts, gradient=True)

            nextdens = density + gradient.T * dx

            zind = np.where((nextdens > self.surfaceCutoff).any(axis=0) & \
                            (density < self.surfaceCutoff))[0]
            zind = np.append(zind, np.where((nextdens < self.surfaceCutoff).any(axis=0) & \
                            (density > self.surfaceCutoff))[0])

            if len(zind) > 0:
                # surfacePoints.extend(pts[zind])
                for i in zind:
                    point = pts[i]

                    if refine:
                        # newton-raphson step in each direction
                        nrstep = (density[i] - self.surfaceCutoff) / gradient[i]
                        # take the one that passes the surface threshold
                        mndx = np.argmin(np.abs(nrstep))
                        # perform step TODO: I though it should have been -= (but doesn't work)
                        point[mndx] += nrstep[mndx]

                    surfacePoints.append(point)

        return np.array(surfacePoints)

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

    def calculate_distance_vector(self, x1, x2):
        dst = x1 - x2
        dst -= self.pbc * np.round(dst / self.pbc)

        return dst

    def calculate_distance(self, x1, x2):
        distance = np.array(self.calculate_distance_vector(x1, x2), ndmin=2)
        return np.linalg.norm(distance, axis=1)

default_zeta = {
    'O' : 2.5, # Angstrom
    'H' : 0.,
    }
