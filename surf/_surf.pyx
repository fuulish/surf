
import cython
import numpy as np
cimport numpy as np
from libc.math cimport round, pow

DTYPE = np.float32
ctypedef np.float32_t DTYPE_float

DTYPE = np.float64
ctypedef np.float64_t DTYPE_double

ITYPE = np.int
ctypedef np.int_t DTYPE_int

cdef extern double c_coarse_grained_density( double *mepos, double *pos, double *zeta, long int natoms,
                                   double *pbc, long int periodic, long int calc_grad, double *grad );

@cython.boundscheck(False)
@cython.wraparound(False)
def coarse_grained_density(
        np.ndarray[double, ndim=1, mode="c"] mepos not None,
        np.ndarray[double, ndim=1, mode="c"] pos not None,
        np.ndarray[double, ndim=1, mode="c"] zeta not None,
        long int natoms,
        np.ndarray[double, ndim=1, mode="c"] pbc not None,
        long int periodic,
        long int calc_grad,
        np.ndarray[double, ndim=1, mode="c"] grad not None,
        ):

    return c_coarse_grained_density(&mepos[0], &pos[0], &zeta[0], natoms, &pbc[0], periodic, calc_grad, &grad[0])

cdef extern double c_opt_distance_to_surface_gsl( double *init_guess, double *mepos, double *pos,
                                      double *zeta, long int natoms, double surfcut, double *pbc,
                                      long int periodic, double *bnds, double xtol, double ctol )

@cython.boundscheck(False)
@cython.wraparound(False)
def opt_distance_to_surface_gsl(
        np.ndarray[double, ndim=1, mode="c"] init_guess not None,
        np.ndarray[double, ndim=1, mode="c"] mepos not None,
        np.ndarray[double, ndim=1, mode="c"] pos not None,
        np.ndarray[double, ndim=1, mode="c"] zeta not None,
        long int natoms,
        double surfcut,
        np.ndarray[double, ndim=1, mode="c"] pbc not None,
        long int periodic,
        np.ndarray[double, ndim=1, mode="c"] bnds not None,
        double xtol,
        double ctol
        ):

    return c_opt_distance_to_surface_gsl( &init_guess[0], &mepos[0], &pos[0],
           &zeta[0], natoms, surfcut, &pbc[0], periodic, &bnds[0], xtol, ctol)
