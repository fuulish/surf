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
#ifdef HAVE_EINSPLINE
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

cube_t interpolate_cube_trilinear ( cube_t * original, int factor, int periodic )
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
    int m000, m100, m010, m110, m001, m101, m011, m111;
    real x000_dat, x100_dat, x010_dat, x110_dat, x001_dat, x101_dat, x011_dat, x111_dat;
    real xd, yd, zd;
    real dx[DIM];

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

    get_box_volels_pointer(original, dx);

    count = 0;
    fctcnt = 0;

    int xup, yup, zup;

    for ( i=0; i<fine.n[0]; i++)
        for ( j=0; j<fine.n[1]; j++)
            for ( k=0; k<fine.n[2]; k++)
            {
                if ( fctcnt == factor )
                    fctcnt = 0;

                // check here, should this be roundf?
                if ( fctcnt == 0 ) {
                    x = (int) floor((float)i/rfct);
                    y = (int) floor((float)j/rfct);
                    z = (int) floor((float)k/rfct);
                }

                index = z + original->n[2] * ( y + original->n[1] * x);

                xup = x+1;
                yup = y+1;
                zup = z+1;

                if ( periodic ) {
                    if ( zup == original->n[2] )
                        periodify_indices ( &zup, &(original->n[2]), &zup, 1);

                    if ( yup == original->n[1] )
                        periodify_indices ( &yup, &(original->n[1]), &yup, 1);

                    if ( xup == original->n[0] )
                        periodify_indices ( &xup, &(original->n[0]), &xup, 1);

                }

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

                // this will contain the value of the voxel below for non-periodic calculation
                // and will not contain any reasonable values

                if ( periodic ) {
                    x000_dat = original->voxels[x000].data;
                    x100_dat = original->voxels[x100].data;
                    x010_dat = original->voxels[x010].data;
                    x110_dat = original->voxels[x110].data;
                    x001_dat = original->voxels[x001].data;
                    x101_dat = original->voxels[x101].data;
                    x011_dat = original->voxels[x011].data;
                    x111_dat = original->voxels[x111].data;
                }
                else {
                    if ( ( zup == original->n[2] ) || ( yup == original->n[1] ) || ( xup == original->n[0] ) ) {
                        m100 = z + original->n[2] * ( y + original->n[1] * (x-1) );
                        m010 = z + original->n[2] * ( (y-1) + original->n[1] * x );
                        m110 = z + original->n[2] * ( (y-1) + original->n[1] * (x-1) );
                        m001 = (z-1) + original->n[2] * ( y + original->n[1] * x );
                        m101 = (z-1) + original->n[2] * ( y + original->n[1] * (x-1) );
                        m011 = (z-1) + original->n[2] * ( (y-1) + original->n[1] * x );
                        m111 = (z-1) + original->n[2] * ( (y-1) + original->n[1] * (x-1) );

                        x000_dat = original->voxels[x000].data;
                        x100_dat = original->voxels[x000].data;
                        x010_dat = original->voxels[x000].data;
                        x110_dat = original->voxels[x000].data;
                        x001_dat = original->voxels[x000].data;
                        x101_dat = original->voxels[x000].data;
                        x011_dat = original->voxels[x000].data;
                        x111_dat = original->voxels[x000].data;

                        //now extrapolate the new data from the gradient to the voxel below the current voxel

                        // if ( ( zup == original->n[2] ) && ( yup == original->n[1] ) && ( xup == original->n[0] ) )
                            // no voxel data needs to be changed, because everything needs to be extrapolated

                        if ( ( zup == original->n[2] ) && ( yup == original->n[1] ) && ( xup != original->n[0] ) )
                            x100_dat = original->voxels[x100].data;

                        else if ( ( zup == original->n[2] ) && ( yup != original->n[1] ) && ( xup == original->n[0] ) )
                            x010_dat = original->voxels[x010].data;

                        else if ( ( zup != original->n[2] ) && ( yup == original->n[1] ) && ( xup == original->n[0] ) )
                            x001_dat = original->voxels[x001].data;

                        else if ( ( zup != original->n[2] ) && ( yup != original->n[1] ) && ( xup == original->n[0] ) )
                            x011_dat = original->voxels[x011].data;

                        else if ( ( zup != original->n[2] ) && ( yup == original->n[1] ) && ( xup != original->n[0] ) )
                            x101_dat = original->voxels[x000].data;

                        else if ( ( zup == original->n[2] ) && ( yup != original->n[1] ) && ( xup != original->n[0] ) )
                            x110_dat = original->voxels[x110].data;

                        else if ( ( zup != original->n[2] ) && ( yup != original->n[1] ) && ( xup != original->n[0] ) )
                            x111_dat = original->voxels[x000].data;

                        if ( zup == original->n[2] ) {
                            real delz = original->voxels[x000].data - original->voxels[m001].data;

                            x001_dat += 2.*delz;
                            x101_dat += 2.*delz;
                            x011_dat += 2.*delz;
                            x111_dat += 2.*delz;

                        }

                        if ( yup == original->n[1] ) {
                            real dely = original->voxels[x000].data - original->voxels[m010].data;

                            x010_dat += 2.*dely;
                            x110_dat += 2.*dely;
                            x011_dat += 2.*dely;
                            x111_dat += 2.*dely;
                        }

                        if ( xup == original->n[0] ) {
                            real delx = original->voxels[x000].data - original->voxels[m100].data;

                            x100_dat += 2.*delx;
                            x110_dat += 2.*delx;
                            x101_dat += 2.*delx;
                            x111_dat += 2.*delx;
                        }
                    }
                    else {
                        x000_dat = original->voxels[x000].data;
                        x100_dat = original->voxels[x100].data;
                        x010_dat = original->voxels[x010].data;
                        x110_dat = original->voxels[x110].data;
                        x001_dat = original->voxels[x001].data;
                        x101_dat = original->voxels[x101].data;
                        x011_dat = original->voxels[x011].data;
                        x111_dat = original->voxels[x111].data;
                    }
                }

                c00 = x000_dat * ( 1 - xd ) + x100_dat * xd;
                c10 = x010_dat * ( 1 - xd ) + x110_dat * xd;
                c01 = x001_dat * ( 1 - xd ) + x101_dat * xd;
                c11 = x011_dat * ( 1 - xd ) + x111_dat * xd;

                c0 = c00 * ( 1 - yd) + c10 * yd;
                c1 = c01 * ( 1 - yd) + c11 * yd;

                c = c0 * ( 1 - zd ) + c1 * zd;

                fine.voxels[count].data = c;

                count++;
                fctcnt ++;
            }

    return fine;
}

#ifdef HAVE_EINSPLINE
cube_t interpolate_cube_bsplines ( cube_t * original, int factor, int periodic )
{
    int i, j, k, indx;

    real cboxv[DIM][DIM];
    float x, y, z, val;
    int cn[DIM];
    cube_t fine;
    real fdx[DIM];

    real rfct = (real) factor;

    // create 3D spline object
	UBspline_3d_s *spline_3d_xyz = get_cube_bsplines ( original, periodic );

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

    destroy_Bspline(spline_3d_xyz);

    return fine;
}

UBspline_3d_s * get_cube_bsplines ( cube_t * cube, int periodic )
{
    int i, j, k;
    real dx[DIM];
    get_box_volels_pointer(cube, dx);
    Ugrid x_grid, y_grid, z_grid;

    // set x_grid, y_grid, z_grid to correct values
    // other two data entries (delta, delta_inv) are considered private

    x_grid.start = cube->origin[0];
    x_grid.end = cube->origin[0]+cube->n[0]*dx[0];
    x_grid.num = cube->n[0];

    y_grid.start = cube->origin[1];
    y_grid.end = cube->origin[1]+cube->n[1]*dx[1];
    y_grid.num = cube->n[1];

    z_grid.start = cube->origin[2];
    z_grid.end = cube->origin[2]+cube->n[2]*dx[2];
    z_grid.num = cube->n[2];

    BCtype_s xBC, yBC, zBC;

    // boundary conditions (periodic)
    if ( periodic ) {
	    xBC = set_boundary_conditions_bsplines ( PERIODIC, PERIODIC, cube->origin[0], cube->origin[0] );
	    yBC = set_boundary_conditions_bsplines ( PERIODIC, PERIODIC, cube->origin[1], cube->origin[1] );
	    zBC = set_boundary_conditions_bsplines ( PERIODIC, PERIODIC, cube->origin[2], cube->origin[2] );
    }
    else {
        // check here, need to figure out correct way of getting non-periodic BSpline, this looks okayish so far
	    xBC = set_boundary_conditions_bsplines ( NATURAL, NATURAL, cube->origin[0], cube->origin[0] );
	    yBC = set_boundary_conditions_bsplines ( NATURAL, NATURAL, cube->origin[1], cube->origin[1] );
	    zBC = set_boundary_conditions_bsplines ( NATURAL, NATURAL, cube->origin[2], cube->origin[2] );
    }

    float * data = (float *) malloc(cube->nvoxels * sizeof ( float ));

    // reassign data to new array that EINSPLINE can work with
    // check here, in general, maybe we should also use just a contiguous array of data points instead of array of structure voxel_t
    for ( i=0; i<cube->nvoxels; i++ )
        data[i] = cube->voxels[i].data;

    // create 3D spline object
	UBspline_3d_s *spline_3d_xyz = create_UBspline_3d_s(x_grid, y_grid, z_grid, xBC, yBC, zBC, data);

    free ( data );

    return spline_3d_xyz;
}

BCtype_s set_boundary_conditions_bsplines ( bc_code lp, bc_code rp, float lVal, float rVal )
{
    BCtype_s bc;

    bc.lCode = lp;
    bc.rCode = rp;

    bc.lVal = lVal;
    bc.rVal = rVal;

    return bc;
}

#endif

cube_t local_interpolation ( cube_t *cube, real *point, int lint, int interpolkind, int ninterpol, char *outputprefix, real *pbc, int periodic )
{
    // find which voxel that point belongs to

    int ix[DIM];
    int ivx;
    real dx[DIM];

    get_box_volels_pointer ( cube, dx );

    get_index_triple ( ix, point, pbc, cube->origin, cube->n, dx, periodic );
    ivx = get_index ( cube->n, ix[0], ix[1], ix[2] );

    // create fake box around that point
    real origin[DIM];
    int cn[DIM];

    int l, m, n;

    for ( l=0; l<DIM; l++ ) {
        origin[l] = cube->voxels[ivx].coords[l] - lint * dx[l];
        cn[l] = ( lint*2 + 1 );
    }

    cube_t cutcube, fine;
    cutcube = initialize_cube ( origin, cube->boxv, cn, cube->atoms, cube->natoms );

    // assign the data from the original cube file into small cutout

    int mnx, mxx, mny, mxy, mnz, mxz;

    // check here, this was just a quick work-around
    mnx = ix[0] - lint;
    mxx = ix[0] + lint + 1;

    mny = ix[1] - lint;
    mxy = ix[1] + lint + 1;

    mnz = ix[2] - lint;
    mxz = ix[2] + lint + 1;

    // check here, maybe just don't use the periodicity at all
    if ( !(periodic) ) {
        if ( ( mnx < 0 ) || ( mny < 0 ) || (mnz < 0 ) || ( mxx >= cube->n[0] ) || ( mxy >= cube->n[1] ) || ( mxz >= cube->n[2] ) ) {
            printf("Would segfault, because non-periodic and local interpolation outside of the cube boundaries, so stopping now!\n");
            exit ( 1 );
        }
    }

    int oindx, nindx;

    int cnl = 0;
    int cnm = 0;
    int cnn = 0;

    int il, im, in;

    cnl = 0;
    for ( l=mnx; l<mxx; l++ ) {
        cnm = 0;
        for ( m=mny; m<mxy; m++ ) {
            cnn = 0;
            for ( n=mnz; n<mxz; n++ ) {
                // check here!!!
                if ( periodic ) {
                    periodify_indices ( &il, &(cube->n[0]), &l, 1 );
                    periodify_indices ( &im, &(cube->n[1]), &m, 1 );
                    periodify_indices ( &in, &(cube->n[2]), &n, 1 );
                }

                oindx = get_index ( cube->n, il, im, in );
                nindx = get_index ( cutcube.n, cnl, cnm, cnn );

                cutcube.voxels[nindx].data = cube->voxels[oindx].data;

                cnn++;
            }
            cnm++;
        }
        cnl++;
    }

#ifdef DEBUG
    char tmp[MAXSTRLEN];
    sprintf(tmp, "%s%i_%s", outputprefix, r, "test_local-non-interpolation.cube");
    write_cubefile(tmp, &cutcube);
#endif
    // interpolate in that region

    // check here, local interpolation is always done non-periodically right now (that is no problem, because the cutout should be assigned according to periodicity and local, periodic interpolation will lead to artifacts)
    if ( interpolkind == INTERPOLATE_TRILINEAR )
        fine = interpolate_cube_trilinear ( &cutcube, ninterpol, 0 );
#ifdef HAVE_EINSPLINE
    else if ( interpolkind == INTERPOLATE_BSPLINES )
        fine = interpolate_cube_bsplines ( &cutcube, ninterpol, 0 );
#endif

#ifdef DEBUG
    sprintf(tmp, "%s%i_%s", outputprefix, r, "test_local-interpolation.cube");
    write_cubefile(tmp, &fine);
#endif

    free ( cutcube.atoms );
    free ( cutcube.voxels );

    return fine;

}
