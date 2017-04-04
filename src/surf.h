/*

SURF is a C program to generate and analyse instantaneous liquid interfaces.
Copyright 2015 Frank Uhlig (uhlig.frank@gmail.com)

This file is part of SURF.

SURF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SURF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SURF.  If not, see <http://www.gnu.org/licenses/>.
*/

cube_t instant_surface_periodic ( int * mask, atom_t * atoms, int inpnatoms, real *zeta, real surfcut, int output, char * outputprefix, real * pbc, real resolution, real accuracy, int provide_box, real * origincube, int * ncube, real boxvcube[DIM][DIM], int periodic, int provide_mask );
real ** get_2d_representation_ils ( int * nsurf, int ** drctn, real ** grad, cube_t * surface, real surfcut, int newsurf, int * surf_inds, int direction, real * area, int periodic );
real get_distance_to_surface ( int * mnnd, int nsurf, real ** surfpts, int * direction, real * grad, atom_t * atoms, int * refmask, int nref, int natoms, real * pbc, int output, char * opref, real surfcut, int periodic );
int check_if_surface_voxel ( int * upper, int * lower, real * tmpdt, cube_t * surface, int * ix, int direction, real surfcut, int periodic );
real get_bulk_volume ( cube_t * surface, real surfcut );
double get_coarse_grained_density( double *mepos, int * mask, atom_t * atoms, real *zeta, real * pbc, int periodic, double * grad );
#ifdef HAVE_NLOPT
real get_opt_distance_to_surface( real *init_guess, real *mepos, int *mask, atom_t * atoms, real *zeta, real surfcut, real *pbc, int periodic, real *bnds, real xtol, real ctol );
#endif
