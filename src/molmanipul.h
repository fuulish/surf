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

void distance_matrix(atom_t * atoms, int * mask, int natoms, real ** distmat, int periodic, real *pbc);
int ** guess_fragments(atom_t * atoms, int natoms, int * mask, int * numfrags, real blfudge, int periodic, real *pbc);
void recurse_concatenate ( char * frags, int ** bondinds, int * countbonds, int atomind, int * isprocessed, int * natomsfrag);
void get_center_of_mass ( real * com, atom_t * atoms, int * mask, int natoms );
void move_atoms ( atom_t * atoms, int * mask, int natoms, real * mover );
