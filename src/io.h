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

#include <stdio.h>
#include <xdrfile/xdrfile.h>

void parse_input_file(input_t *, char *);
void set_input_value(input_t *, char *, char *);
void read_cubefile(char *, cube_t *, int );
void read_cubefile_data(char *, real * );
void write_cubefile(char *, cube_t *);
void write_cubefile_array(char * filename, cube_t * cube, real data[]);
void write_cubefile_indices(char *, cube_t *, int *); // int);
void write_cubefile_indices_plain(char *, cube_t *, int *);
void write_cubefile_offset(char *, cube_t *, int);
void allocate_cube(voxel_t **, int);
void allocate_atoms(atom_t **, int);
void set_input_defaults(input_t *);
int * get_mask(char *, char *, int, atom_t * atoms, int natoms);
int read_xmol(char *, atom_t **);
int check_input_dependencies(input_t *);
int replicate_structure(cube_t * cube, real pbc[], int nreplicas[]);
int add_atom(cube_t * cube, char * kind, char * symdex, real * center);
void write_cubefile_indices_number(char * filename, cube_t * cube, int * indices, int number);
real integrate_cube_number(int number, int * indices, cube_t * cube);
void write_chgcar(char * filename, cube_t * cube);
void parse_cmdline(input_t * inppar, char ** argv, int argc);
void display_usage(void);
void xmolreader(FILE * fxmol, int bytelen, int snap, atom_t * atoms, int natoms);
int xmol_snap_bytesize(FILE * fxmol);
void read_xtr_forward (XDRFILE * xd_read, int frwrd, atom_t * atoms, int natoms);
void write_xyz ( FILE * xyz, atom_t * atoms, int natoms, char * comment );
void read_roots_file ( char * rootsfile, int * nx, int * ny );
void get_dimensions_roots ( char * fname, int * n1, int * n2 );
void write_matrix_real_2d_to_file ( char * fname, real ** mat, int n1, int n2 );
void write_matrix_real_2d_to_file_w_cont_spacing ( char * fname, real ** mat, int n1, int n2, real *origin, real d1, real d2 );
int assign_task ( char * task );
