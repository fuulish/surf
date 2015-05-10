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

#include "declarations.h"
#include "types.h"
#include "constants.h"
#include "io.h"
#include "utils.h"
#include "atom_param.h"
#include "cube.h"
#include "molmanipul.h"
#include "errors.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef EINSPLINE
#include <einspline/bspline.h>
#endif

int tstart;
int tstop;

cube_t initialize_cube(real origin[DIM], real boxv[DIM][DIM], int n[DIM], atom_t * atoms, int natoms)
{
    int i, j, k, l;
    cube_t cube;

    strcpy(cube.title, "interpolated cube\n");
    strcpy(cube.comment, "no comment\n");

    for ( i=0; i<DIM; i++)
        cube.origin[i] = origin[i];

    for ( i=0; i<DIM; i++ )
        for ( j=0; j<DIM; j++ )
            cube.boxv[i][j] = boxv[i][j];

    cube.dv = boxv[0][0] * ( boxv[1][1] * boxv[2][2] - boxv[2][1] * boxv[1][2] ) 
             - boxv[0][1] * ( boxv[1][0] * boxv[2][2] - boxv[2][0] * boxv[1][2] )
             + boxv[0][2] * ( boxv[1][0] * boxv[2][1] - boxv[2][0] * boxv[1][1] );

    if ( (sqrt(boxv[0][0] * boxv[1][0] + boxv[0][1] * boxv[1][1] + boxv[0][2] * boxv[1][2]) == 0) && 
         (sqrt(boxv[0][0] * boxv[2][0] + boxv[0][1] * boxv[2][1] + boxv[0][2] * boxv[2][2]) == 0) && 
         (sqrt(boxv[1][0] * boxv[2][0] + boxv[1][1] * boxv[2][1] + boxv[1][2] * boxv[2][2]) == 0))
    {
        cube.orthogonal = 1;
    }
    else
    {
        cube.orthogonal = 0;
    }

    if ( (sqrt(sqr(boxv[0][0]) + sqr(boxv[0][1]) + sqr(boxv[0][2])) - sqrt(sqr(boxv[1][0]) + sqr(boxv[1][1]) + sqr(boxv[1][2])) < EPS) &&
         (sqrt(sqr(boxv[0][0]) + sqr(boxv[0][1]) + sqr(boxv[0][2])) - sqrt(sqr(boxv[2][0]) + sqr(boxv[2][1]) + sqr(boxv[2][2])) < EPS))
    {
        cube.cubic = 1;
    }
    else
    {
        cube.cubic = 0;
    }

    cube.nvoxels = 1;

    for (i=0; i<DIM; i++)
    {
        cube.n[i] = n[i];
        cube.nvoxels *= n[i];
    }

    cube.natoms = natoms;
    cube.atoms = (atom_t*) malloc(cube.natoms * (sizeof(atom_t)));
    if ( NULL == ( cube.voxels = (voxel_t*) malloc(cube.nvoxels * (sizeof(voxel_t))))) {
        print_error ( OUT_OF_MEMORY, "allocation of cube voxels" );
        exit ( OUT_OF_MEMORY );
    }

    for (i=0; i<cube.natoms; i++)
    {
        char buffer[MAXSTRLEN];
        sprintf(buffer, "%i", atoms[i].number);

        assign_atom_parameters("index", buffer, &(cube.atoms[i]));
        
        for (j=0; j<DIM; j++)
            cube.atoms[i].coords[j] = atoms[i].coords[j];
    }

    int cuben[DIM];
    for (k=0; k<DIM; k++)
        cuben[k] = 0;

    for ( i=0; i<cube.nvoxels; i++ )
    {
        cube.voxels[i].data = ZERO;

        for (k=0; k<DIM; k++)
        {
            cube.voxels[i].coords[k] = ZERO;
            for (l=0; l<DIM; l++)
                cube.voxels[i].coords[l] += cube.boxv[k][l] * cuben[l];

            cube.voxels[i].coords[k] += cube.origin[k];
        }

        cuben[DIM-1]++;

        switch (DIM)
        {
            case 3:
                {
                    if (cuben[2] == cube.n[2])
                    {
                        cuben[2] = 0;
                        cuben[1]++;
                    }

                    if (cuben[1] == cube.n[1])
                    {
                        cuben[1] = 0;
                        cuben[0]++;
                    }
                    break;
                };
            case 2:
                {
                    if (cuben[1] == cube.n[1])
                    {
                        cuben[1] = 0;
                        cuben[0]++;
                    }
                    break;
                }
        }
    }

    return cube;
}

int * cubes_larger_than(real cutoff, cube_t * cube)
{
    int i, k;
    int * indices;
    
    indices = (int *) malloc((cube->nvoxels+1) * sizeof(int));

    k = 0;

    for (i=0; i<cube->nvoxels; i++)
    {
        if (cube->voxels[i].data > cutoff)
        {
            indices[k] = i;
            k++;
        }
    }

    indices[k] = -1;

    return indices;
}

int * invert_indices(int nvoxels, int * indices)
{
    int i, j, k;
    int * inverted;

    inverted = (int *) malloc((nvoxels+1) * sizeof(int));

    k = 0;
    j = 0;

    // this only works if the indices array is sorted in ascending order
    for (i=0; i<nvoxels; i++)
    {
        if ( i == indices[k] )
        {
            k++;
        }
        else
        {
            inverted[j] = i;
            j++;
        }
    }

    inverted[j] = -1;

    return inverted;
}

void get_cell_pointer(cube_t * cube, real * cell)
{
    int i, j;
    real dx[DIM];
    
    for ( i=0; i<DIM; i++)
        dx[i] = ZERO;

    for ( i=0; i<DIM; i++ )
    {
        for (j=0; j<DIM; j++)
        {
            dx[i] += cube->boxv[i][j];
        }

        cell[i] = cube->n[i] * dx[i];
    }

}

void get_box_volels_pointer(cube_t * cube, real * dx)
{
    int i, j;
    
    for ( i=0; i<DIM; i++)
        dx[i] = ZERO;

    for ( i=0; i<DIM; i++ )
        for (j=0; j<DIM; j++)
            dx[i] += cube->boxv[i][j];

}

real * get_box_volels(cube_t * cube)
{
    int i, j;
    real * dx;
    
    dx = (real *) malloc(DIM * sizeof(real));

    for ( i=0; i<DIM; i++)
        dx[i] = ZERO;

    for ( i=0; i<DIM; i++ )
        for (j=0; j<DIM; j++)
            dx[i] += cube->boxv[i][j];

    return dx;

}

void get_box_areas_pointer (real * da, cube_t * cube, real * dx )
{
    int i, j;

    for ( i=0; i<DIM; i++ )
        da[i] = ONE;

    for ( i=0; i<DIM; i++ ) {
        for ( j=0; j<DIM; j++ ) {
            da[i] *= dx[j];

            if ( i == j )
                da[i] /= dx[j];
        }
    }

}

cube_t interpolate_cube_trilinear ( cube_t * original, int factor )
{
    int i, j, k;
    cube_t fine;
    real cboxv[DIM][DIM];
    int cn[DIM];
    int count, valcounter;
    int fctcnt;
    int index;
    int x, y, z;
    int x000, x100, x010, x110, x001, x101, x011, x111;
    real xd, yd, zd;

    real rfct = (real) factor;

    real c00, c10, c01, c11, c0, c1, c;

    for ( i=0; i<DIM; i++ )
    {
        cn[i] = original->n[i]*factor;
        for ( j=0; j<DIM; j++ )
        {
            cboxv[i][j] = original->boxv[i][j] / rfct;
        }
    }

    fine = initialize_cube(original->origin, cboxv, cn, original->atoms, original->natoms);

    count = 0;
    fctcnt = 0;

    int xup, yup, zup;

    for ( i=0; i<fine.n[0]; i++)
        for ( j=0; j<fine.n[1]; j++)
            for ( k=0; k<fine.n[2]; k++)
            {
                if ( fctcnt == factor )
                    fctcnt = 0;

                if ( fctcnt == 0 ) {
                    x = (int) floor((float)i/rfct);
                    y = (int) floor((float)j/rfct);
                    z = (int) floor((float)k/rfct);
                }

                index = z + original->n[2] * ( y + original->n[1] * x);

                xup = x+1;
                yup = y+1;
                zup = z+1;

                if ( zup == original->n[2] )
                    periodify_indices ( &zup, &(original->n[2]), &zup, 1);

                if ( yup == original->n[1] )
                    periodify_indices ( &yup, &(original->n[1]), &yup, 1);

                if ( xup == original->n[0] )
                    periodify_indices ( &xup, &(original->n[0]), &xup, 1);

                xd = (fine.voxels[count].coords[0] - original->voxels[index].coords[0]) / original->boxv[0][0];
                yd = (fine.voxels[count].coords[1] - original->voxels[index].coords[1]) / original->boxv[1][1];
                zd = (fine.voxels[count].coords[2] - original->voxels[index].coords[2]) / original->boxv[2][2];

                x000 = z + original->n[2] * ( y + original->n[1] * x );
                x100 = z + original->n[2] * ( y + original->n[1] * (xup) );
                x010 = z + original->n[2] * ( (yup) + original->n[1] * x );
                x110 = z + original->n[2] * ( (yup) + original->n[1] * (xup) );
                x001 = (zup) + original->n[2] * ( y + original->n[1] * x );
                x101 = (zup) + original->n[2] * ( y + original->n[1] * (xup) );
                x011 = (zup) + original->n[2] * ( (yup) + original->n[1] * x );
                x111 = (zup) + original->n[2] * ( (yup) + original->n[1] * (xup) );

                c00 = original->voxels[x000].data * ( 1 - xd ) + original->voxels[x100].data * xd;
                c10 = original->voxels[x010].data * ( 1 - xd ) + original->voxels[x110].data * xd;
                c01 = original->voxels[x001].data * ( 1 - xd ) + original->voxels[x101].data * xd;
                c11 = original->voxels[x011].data * ( 1 - xd ) + original->voxels[x111].data * xd;

                c0 = c00 * ( 1 - yd) + c10 * yd;
                c1 = c01 * ( 1 - yd) + c11 * yd;

                c = c0 * ( 1 - zd ) + c1 * zd;

                fine.voxels[count].data = c;

                count++;
                fctcnt ++;
            }

    return fine;
}

#ifdef EINSPLINE
cube_t interpolate_cube_bsplines ( cube_t * original, int factor )
{

    int i, j, k, indx, einindx;

    real cboxv[DIM][DIM];
    float x, y, z, val;
    int cn[DIM];
    cube_t fine;
    real dx[DIM], fdx[DIM];

    get_box_volels_pointer(original, dx);
    real rfct = (real) factor;

    Ugrid x_grid, y_grid, z_grid;

    // set x_grid, y_grid, z_grid to correct values
    // other two data entries (delta, delta_inv) are considered private

    x_grid.start = original->origin[0];
    x_grid.end = original->origin[0]+original->n[0]*dx[0];
    x_grid.num = original->n[0];

    y_grid.start = original->origin[1];
    y_grid.end = original->origin[1]+original->n[1]*dx[1];
    y_grid.num = original->n[1];

    z_grid.start = original->origin[2];
    z_grid.end = original->origin[2]+original->n[2]*dx[2];
    z_grid.num = original->n[2];

    // boundary conditions (periodic)
    // what is the meaning of the last two entries? left and right boundary conditions values (not sure what exactly)
    // maybe we should set boundary values to original->origin[i], but this is just a guess
	BCtype_s xBC = {PERIODIC, PERIODIC , original->origin[0], original->origin[0]};
	BCtype_s yBC = {PERIODIC, PERIODIC , original->origin[1], original->origin[1]};
	BCtype_s zBC = {PERIODIC, PERIODIC , original->origin[2], original->origin[2]};

    float * data = (float *) malloc(original->nvoxels * sizeof ( float ));

    // reassign data to new array that EINSPLINE can work with
    // check here, in general, maybe we should also use just a contiguous array of data points instead of array of structure voxel_t
    for ( i=0; i<original->nvoxels; i++ )
        data[i] = original->voxels[i].data;

    // create 3D spline object
	UBspline_3d_s *spline_3d_xyz = create_UBspline_3d_s(x_grid, y_grid, z_grid, xBC, yBC, zBC, data);

    // initialize new and fine cube (first, set dimensions, second initialize cube)
    for ( i=0; i<DIM; i++ ) {
        cn[i] = original->n[i]*factor;
        for ( j=0; j<DIM; j++ )
            cboxv[i][j] = original->boxv[i][j] / rfct;
    }

    fine = initialize_cube(original->origin, cboxv, cn, original->atoms, original->natoms);

    get_box_volels_pointer(&fine, fdx);

    // now evaluate spline at new points

    for ( i=0; i<fine.n[0]; i++ )
        for ( j=0; j<fine.n[1]; j++ )
            for ( k=0; k<fine.n[2]; k++ ) {

                x = fine.origin[0] + i* fdx[0];
                y = fine.origin[1] + j* fdx[1];
                z = fine.origin[2] + k* fdx[2];

                val = ZERO;
                eval_UBspline_3d_s ( spline_3d_xyz, x, y, z, &val );

                indx = get_index ( fine.n, i, j, k );
                fine.voxels[indx].data = val;
            }

    free ( data );

    destroy_Bspline(spline_3d_xyz);

    return fine;
}
#endif
