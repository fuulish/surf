
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
