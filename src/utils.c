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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "utils.h"

void put_zeros_real_1d ( real * data, int len )
{
    int i;

    for ( i=0; i<len; i++ )
        data[i] = ZERO;
}

void put_zeros_real_2d ( real ** data, int len1, int len2 )
{
    int i, j;

    for ( i=0; i<len1; i++ )
        for ( j=0; j<len2; j++ )
            data[i][j] = ZERO;
}

void put_zeros_int_1d ( int * data, int len )
{
    int i;

    for ( i=0; i<len; i++ )
        data[i] = ZERO;
}

void put_zeros_int_2d ( int ** data, int len1, int len2 )
{
    int i, j;

    for ( i=0; i<len1; i++ )
        for ( j=0; j<len2; j++ )
            data[i][j] = ZERO;
}

// real * allocate_matrix_real_1d ( int dim )
// {
//     int i;
//     real * mat;
//
//     mat = (real *) malloc (dim * sizeof(real));
//
//     return mat;
// }
//
// void free_matrix_real_1d ( real * mat )
// {
//     int i;
//
//     free ( mat );
// }

real ** allocate_matrix_real_2d ( int dim1, int dim2 )
{
    int i;
    real ** mat;

    mat = (real **) malloc (dim1 * sizeof(real *));

    for ( i=0; i<dim1; i++ )
        mat[i] = (real *) calloc( dim2, sizeof(real));

    return mat;
}

void free_matrix_real_2d ( real ** mat, int dim1 )
{
    int i;
    for ( i=0; i<dim1; i++ )
        free(mat[i]);

    free(mat);
}

int ** allocate_matrix_int_2d ( int dim1, int dim2 )
{
    int ** mat;
    int i;

    mat = (int **) malloc (dim1 * sizeof(int *));

    for ( i=0; i<dim1; i++ )
        mat[i] = (int *) calloc( dim2, sizeof(int));

    return mat;
}
void free_matrix_int_2d ( int ** mat, int dim1 )
{
    int i;
    for ( i=0; i<dim1; i++ )
        free(mat[i]);

    free(mat);
}

// void qsort_real_indices (

int find_minimum_d(real * values, int * mask)
{
    int i, minindex;
    real minimum;

    minimum = values[0];
    minindex = 0;

    i = 0;
    while ( mask[i] != -1 )
    {
        if (values[i] < minimum)
        {
            minimum = values[i];
            minindex = i;
        }

        i++;
    }

    return minindex;
}

real find_minimum_1d_real ( int * mini, real * values, int maxi )
{
    int i;
    real min;

    min = values[0];
    *mini = 0;

    /* check here, there is a more elegant way of doing it, sth like ( condition ) || ( assignment ) */
    for ( i=0; i<maxi; i++ ) {
        if ( values[i] < min ) {
            min = values[i];
            *mini = i;
        }
    }

    return min;
}

real find_maximum_1d_real ( int * mini, real * values, int maxi )
{
    int i;
    real max;

    max = values[0];
    *mini = 0;

    for ( i=0; i<maxi; i++ ) {
        if ( values[i] > max ) {
            max = values[i];
            *mini = i;
        }
    }

    return max;
}

real find_minimum_2d_real ( int * mini, int * minj, real ** values, int maxi, int maxj )
{
    int i, j;
    real min;

    min = values[0][0];
    *mini = 0;
    *minj = 0;

    for ( i=0; i<maxi; i++ )
        for ( j=0; j<maxj; j++ )
            if ( values[i][j] < min ) {
                min = values[i][j];

                *mini = i;
                *minj = j;
            }

    return min;
    // printf("%i %i\n", *mini, *minj);
}

real find_maximum_2d_real ( int * mini, int * minj, real ** values, int maxi, int maxj )
{
    int i, j;
    real max;

    max = values[0][0];
    *mini = 0;
    *minj = 0;

    for ( i=0; i<maxi; i++ )
        for ( j=0; j<maxj; j++ )
            if ( values[i][j] > max ) {
                max = values[i][j];

                *mini = i;
                *minj = j;
            }

    return max;
    // printf("%i %i\n", *mini, *minj);
}

real get_distance_periodic (real * coord1, real * coord2, real * pbc)
{
    int i;
    real distance = ZERO;
    real rxij;

    for ( i=0; i<DIM; i++ ) {
        rxij = coord1[i] - coord2[i];
        rxij -= pbc[i] * roundf( rxij / pbc[i] );

        distance += sqr ( rxij );
    }

    distance = sqrt ( distance );

    return distance;
}

real get_distance_vector_periodic (real * dist, real * coord1, real * coord2, real * pbc)
{
    int i;
    real distance = ZERO;
    real rxij;

    for ( i=0; i<DIM; i++ ) {
        rxij = coord1[i] - coord2[i];
        rxij -= pbc[i] * roundf( rxij / pbc[i] );

        dist[i] = rxij;
    }
}

real get_distance_periodic_1d (real coord1, real coord2, real pbc)
{
    real distance = ZERO;
    real rxij;

    rxij = coord1 - coord2;
    rxij -= pbc * roundf( rxij / pbc );

    if ( rxij < 0 )
        distance = -1 * rxij;
    else
        distance = rxij;

    return distance;
}

real lerp ( real v0, real v1, real t )
{
    // interpolates between v0 and v1 according to the parameter t \element [0,1]
    return v0+(v1-v0)*t;
}

real lerp_to_t ( real v0, real v1, real vt )
{
    // interpolates between v0 and v1 according to find parameter t \element [0,1] for given v(t) (vt)
    return ( vt - v0 ) / ( v1 - v0 );
}

real kahan_sum_real (real * data, int len)
{
    int i;

    real sum = ZERO, c = ZERO;
    real y, t;

    for ( i=0; i<len; i++ ) {
        y = data[i] - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

void cub2arr_add ( real * array, cube_t * cube )
{
    int k;

    for ( k=0; k<cube->nvoxels; k++ )
        array[k] += cube->voxels[k].data;

}

void arr2cub_add ( cube_t * cube, real * array )
{
    int k;

    for ( k=0; k<cube->nvoxels; k++ )
        cube->voxels[k].data += array[k];

}

void cub2arr_sub ( real * array, cube_t * cube )
{
    int k;

    for ( k=0; k<cube->nvoxels; k++ )
        array[k] -= cube->voxels[k].data;

}

void arr2cub_sub ( cube_t * cube, real * array )
{
    int k;

    for ( k=0; k<cube->nvoxels; k++ )
        cube->voxels[k].data -= array[k];

}

void cub2arr_mul ( real * array, cube_t * cube )
{
    int k;

    for ( k=0; k<cube->nvoxels; k++ )
        array[k] *= cube->voxels[k].data;

}

void arr2cub_mul ( cube_t * cube, real * array )
{
    int k;

    for ( k=0; k<cube->nvoxels; k++ )
        cube->voxels[k].data *= array[k];

}

void cub2arr_div ( real * array, cube_t * cube )
{
    int k;

    for ( k=0; k<cube->nvoxels; k++ )
        array[k] /= cube->voxels[k].data;

}

void arr2cub_div ( cube_t * cube, real * array )
{
    int k;

    for ( k=0; k<cube->nvoxels; k++ )
        cube->voxels[k].data /= array[k];

}

void arr2arr_cpy ( real * arrout, real * arrinp, int len )
{
    int k;

    for ( k=0; k<len; k++ )
        arrout[k] = arrinp[k];

}

void arr2arr_add ( real * arrout, real * arrinp, int len )
{
    int k;

    for ( k=0; k<len; k++ )
        arrout[k] += arrinp[k];

}

void arr2arr_sub ( real * arrout, real * arrinp, int len )
{
    int k;

    for ( k=0; k<len; k++ )
        arrout[k] -= arrinp[k];

}

void arr2arr_mul ( real * arrout, real * arrinp, int len )
{
    int k;

    for ( k=0; k<len; k++ )
        arrout[k] *= arrinp[k];

}

void arr2arr_div ( real * arrout, real * arrinp, int len )
{
    int k;

    for ( k=0; k<len; k++ )
        arrout[k] /= arrinp[k];

}

real get_distance(real * coord1, real * coord2)
{
    int i;
    real distance = ZERO;

    for ( i=0; i<DIM; i++ )
        distance += sqr(coord1[i] - coord2[i]);

    distance = sqrt(distance);

    return distance;
}

int get_index ( int * nvox, int i, int j, int k )
{
    return ( k + nvox[2] * ( j + nvox[1] * i ) );
}

void bubble_sort_ints(int indices[], int nindices)
{
    int i, j;
    int h;

    for ( j=nindices-1; j>0; j=j-1 )
    {
        for ( i=0; i<j; i = i + 1 )
        {
            if (indices[i] > indices[i+1])
            {
                h = indices[i];
                indices[i] = indices[i+1];
                indices[i+1] = h;
            }
        }
    }

    // check here, why don't we use pointers?
    return;
}

void cross_product_3d ( real * c, real * a, real * b )
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

real dot_product_nd ( real * a, real * b, int len )
{
    int i;
    real c = ZERO;

    for ( i=0; i<len; i++ )
        c += a[i]*b[i];

    return c;
}

void get_index_triple ( int *i, real* coords, real *pbc, real *resolution, int periodic )
{
    int k;
    int n, u;

    for ( k=0; k<DIM; k++ ) {
        i[k] = (int) floor ( coords[k] / resolution[k] );

        if ( i[k] < 0 ) {
            n = pbc[k] / resolution[k];
            i[k] = n - i[k];
        }
    }
}

void periodify_indices ( int * out, int * nvx, int * in, int len)
{
    int k;

    for ( k=0; k<len; k++ )
        if ( in[k] < 0 )
            out[k] = in[k] + nvx[0];
        else if ( in[k] >= nvx[0] )
            out[k] = in[k] - nvx[0];
        else
            out[k] = in[k];
}
