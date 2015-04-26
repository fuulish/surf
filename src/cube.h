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

cube_t initialize_cube(real origin[DIM], real boxv[DIM][DIM], int n[DIM], atom_t *, int);
int * cubes_larger_than(real cutoff, cube_t * cube);
int * invert_indices(int nvoxels, int * indices);
void get_cell_pointer(cube_t * cube, real * cell);
real * get_box_volels(cube_t * cube);
void get_box_volels_pointer(cube_t * cube, real * dx);
