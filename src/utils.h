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

void put_zeros_real_1d ( real * data, int len );
void put_zeros_real_2d ( real ** data, int len1, int len2 );
void put_zeros_int_1d ( int * data, int len );
void put_zeros_int_2d ( int ** data, int len1, int len2 );
real ** allocate_matrix_real_2d ( int dim1, int dim2 );
void free_matrix_real_2d ( real ** mat, int dim1 );
int ** allocate_matrix_int_2d ( int dim1, int dim2 );
void free_matrix_int_2d ( int ** mat, int dim1 );
int find_minimum_d(real *, int *);
real get_distance_periodic (real * coord1, real * coord2, real * pbc);
real get_distance_vector (real * dist, real * coord1, real * coord2);
void get_distance_vector_periodic (real * dist, real * coord1, real * coord2, real * pbc);
real get_distance_periodic_1d (real coord1, real coord2, real pbc);
real get_distance(real * coord1, real * coord2);
real find_minimum_1d_real ( int * mini, real * values, int maxi );
real find_maximum_1d_real ( int * mini, real * values, int maxi );
real find_minimum_2d_real ( int * mini, int * minj, real ** values, int maxi, int maxj );
real find_maximum_2d_real ( int * mini, int * minj, real ** values, int maxi, int maxj );
int get_index ( int * nvox, int i, int j, int k );
real lerp ( real v0, real v1, real t );
real lerp_to_t ( real v0, real v1, real vt );
real kahan_sum_real (real * data, int len);
double kahan_sum_double (double * data, int len);
real kahan_sum_voxels (voxel_t * data, int len);
int get_func_id ( char * func_name );
void cub2arr_cpy ( real * array, cube_t * cube );
void arr2cub_cpy ( cube_t * cube, real * array );
void cub2arr_add ( real * array, cube_t * cube );
void arr2cub_add ( cube_t * cube, real * array );
void cub2arr_sub ( real * array, cube_t * cube );
void arr2cub_sub ( cube_t * cube, real * array );
void cub2arr_mul ( real * array, cube_t * cube );
void arr2cub_mul ( cube_t * cube, real * array );
void cub2arr_div ( real * array, cube_t * cube );
void arr2cub_div ( cube_t * cube, real * array );

void arr2arr_cpy ( real * arrout, real * arrinp, int len );
void arr2arr_add ( real * arrout, real * arrinp, int len );
void arr2arr_sub ( real * arrout, real * arrinp, int len );
void arr2arr_mul ( real * arrout, real * arrinp, int len );
void arr2arr_div ( real * arrout, real * arrinp, int len );
void bubble_sort_ints(int indices[], int nindices);
void cross_product_3d ( real * c, real * a, real * b );
real dot_product_nd ( real * a, real * b, int len );
void get_index_triple ( int *i, real* coords, real *pbc, real *origin, int *n, real *resolution, int periodic );
void periodify_indices ( int * out, int * nvx, int * in, int len);
