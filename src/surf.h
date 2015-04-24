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

cube_t instant_surface_periodic ( int * mask, atom_t * atoms, int inpnatoms, real zeta, real surfcut, int output, char * outputprefix, real * pbc, real resolution, real accuracy, int provide_box, real * origincube, int * ncube, real boxvcube[DIM][DIM], int periodic, int provide_mask );
void get_2d_representation_ils ( cube_t * surface, real ** surf_2d_up, real ** surf_2d_down, real surfcut, int newsurf, int * surf_up_inds, int * surf_down_inds, int direction );
void get_distance_to_surface ( real * disthi, real * distlo, int * inthi, int * intlo, cube_t * surface, real ** surf_2d_up, real ** surf_2d_down, atom_t * atoms, int * refmask, int nref, int natoms, real * pbc, int output, char * opref, int direction, real surfcut );
