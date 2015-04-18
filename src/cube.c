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
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
    cube.voxels = (voxel_t*) malloc(cube.nvoxels * (sizeof(voxel_t)));

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

