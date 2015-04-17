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

        // check here, i think it's correct, but did not follow it through thoroughly
        for (k=0; k<DIM; k++)
        {
            cube.voxels[i].coords[k] = ZERO;
            // check here, because i think we are treating voxels as grid points
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

cube_t instant_surface_periodic ( int * mask, atom_t * atoms, int inpnatoms, real zeta, real surfcut, int output, char * outputprefix, real * pbc, real resolution, real accuracy, int provide_box, real * origincube, int * ncube, real boxvcube[DIM][DIM], int periodic, int provide_mask )
{
    int i, j, natoms;
    int a;
    int * sindices;
    int * bindices;
    real sqzeta;
    real trplzt;
    cube_t surface;
    real prefactor, dummy;

    real orig[DIM];
    real boxv[DIM][DIM];
    real refc[DIM];
    int n[DIM];

    trplzt = 3*zeta;

    /* need input :
     *      - atomic coordinates of "surface atoms"
     *      - resolution of surface grid
     * only go to 3\sigma and shift the gaussian accordingly
     *      - count how many voxels need to be included in each direction and then reduce computational cost
     *
     */

    /*
     * for non-periodic calculation only add so much space around original structure
     * that we do not have more than 3*zeta from the atoms on the outside
     */

    /* we will create an orthogonal box according to periodic boundary conditions and resolution input */
    if ( provide_box ) {
        for ( i=0; i<DIM; i++ ) {

            orig[i] = origincube[i];
            // check here if there is a nicer reference center, e.g., center of actual cube we are working with
            refc[i] = pbc[i] / 2.;
            n[i] = ncube[i];

            for ( j=0; j<DIM; j++ ) {
                boxv[i][j] = boxvcube[i][j];
            }
        }
    }
    else {
        for ( i=0; i<DIM; i++ ) {

            orig[i] = ZERO;
            refc[i] = pbc[i] / 2.;
            n[i] = pbc[i] / resolution;

            for ( j=0; j<DIM; j++ ) {
                if ( i == j )
                    boxv[i][j] = pbc[i] / n[i];
                else
                    boxv[i][j] = ZERO;
            }
        }
    }

    natoms = 0;
    while ( mask[natoms] != -1 )
        natoms++;

    /* initialize surface cube with original atom positions and new cube dimensions */

    surface = initialize_cube(orig, boxv, n, atoms, inpnatoms);

    /* wrap atoms into central box
     * to this end, center point of box is taken as reference
     */

    /* check here and already include the two down the code already in sqzeta and also try float precision */

    sqzeta = sqr(zeta);
    real mttsqzeta = -2. * sqzeta;
    dummy = 2. * PI * sqzeta;
    prefactor = 1. / dummy / (sqrt(dummy));

    /* get maximum distance in terms of number of voxels, according to accuracy criterion */

    i = 1;

    real maxdist = i * boxv[0][0];

    real cubpbc[DIM];
    real cubedge[DIM*2];
    real dx[DIM];
    int addspace[DIM];

    get_cell_pointer ( &surface, cubpbc );
    get_box_volels_pointer ( &surface, dx );

    for ( i=0; i<DIM; i++ ) {
        cubedge[i] = surface.origin[i] - maxdist;
        cubedge[DIM+i] = surface.origin[i] + cubpbc[i] + maxdist;
        // printf("%f %f\n", cubedge[i], cubedge[DIM+i]);
    }

    /* check here, this will only work if our pbc always starts at (0,0,0) */
    for ( i=0; i<DIM; i++ ) {
        // round down to not overestimate things here
        /* check here, not sure which one is the correct one */

        addspace[i] = lroundf ( ( pbc[i] - cubpbc[i] + surface.origin[i] ) / dx[i] );
    }

    real distance;

    int skip;

    // real rxij[DIM];
    // not sure right now what happens if we pass array as private in openmp loop
    // check here, what happens if we define it inside the shared region (it oughta be private to each thread of the team

    // for non-periodic calculation, sort out atoms that are within range of the surface cube

    int * actatoms;
    int nactive = 0;

    int cnt = 0;

    if ( !(periodic ) ) {
        actatoms = ( int * ) malloc ( ( natoms + 1 ) * sizeof ( int ) );
        /* this should be more accurate than the previous version and has little more overhead */

        for ( a=0; a<natoms; a++ ) {
            // printf("%i %i\n", a, mask[a]);
            skip = 0;

            for ( j=0; j<DIM; j++ ) {
                if ( ( atoms[mask[a]].coords[j] < cubedge[j] ) || ( atoms[mask[a]].coords[j] > cubedge[DIM+j] ) ) {
                    skip = 1;
                    break;
                }
            }

            if ( !( skip ) ) {
                // printf("%s %14.8f %14.8f %14.8f\n", atoms[mask[a]].symbol, BOHR*atoms[mask[a]].coords[0], BOHR*atoms[mask[a]].coords[1], BOHR*atoms[mask[a]].coords[2]);

                // printf("%i\n", mask[a]);
                actatoms[cnt] = mask[a];
                cnt++;
            }
            // for testing purposes
            // else {
            //     printf("%5s %14.8f %14.8f %14.8f\n", atoms[mask[a]].symbol, atoms[mask[a]].coords[0], atoms[mask[a]].coords[1], atoms[mask[a]].coords[2]);
            // }

        }
        actatoms[cnt] = -1;
        nactive = cnt;
    }

#if DEBUG
    printf("%i active atoms for surface generation\n", cnt);
#endif

#ifdef OPTSURF
/* FU| check here, this should become the optimized code
 */

    printf("We are using the optimized, but not debugged version of the code with a triple-zeta of %20.8f!!!\n", trplzt);

    real tmpdst;
    int k;
    real cutshft = prefactor * exp( sqr( trplzt ) / (mttsqzeta));

    if ( periodic ) {
#ifdef OPENMP
#pragma omp parallel for default(none) \
    private(i,j,k,a,distance,skip,tmpdst) shared(surface,atoms,natoms,mask,prefactor,mttsqzeta,pbc,maxdist,cubedge,periodic,trplzt,cutshft) \
        schedule(guided, surface.n[2])
    // schedule(dynamic)
    // it would in principle nice to get chunks in the size of z-slices, so maybe guided(cube.n[2]) would work better? so that we get equal hits and misses for the skip calculation (general volume division / number of threads will prolly work equally well
#endif
        for ( i=0; i<surface.nvoxels; i++)
        {
            for ( a=0; a<natoms; a++ ) {
                distance = ZERO;
                skip = 0;

                for ( k=0; k<DIM; k++ ) {
                    tmpdst = get_distance_periodic_1d ( surface.voxels[i].coords[k], atoms[mask[a]].coords[k], pbc[k] );

                    if ( tmpdst > trplzt )
                        skip = 1;
                    else
                        distance += sqr ( tmpdst );

                    if ( skip ) {
                        // HELLO;
                        continue;
                    }
                }

                if ( skip )
                    continue;

                distance = sqrt ( distance );

                if ( distance > trplzt )
                    continue;

                // still need to shift the value by value of gaussian at 3*sigma

                // distance = get_distance_periodic ( &(surface.voxels[i].coords[0]), &(atoms[mask[a]].coords[0]), pbc );

                surface.voxels[i].data += prefactor * exp( sqr( distance ) / (mttsqzeta));
                surface.voxels[i].data -= cutshft;

                }
        }
    }
    else {
#ifdef OPENMP
#pragma omp parallel for default(none) \
    private(i,j,a,distance,skip) shared(surface,atoms,nactive,actatoms,prefactor,mttsqzeta,pbc,maxdist,cubedge,periodic) \
    // schedule(guided)
#endif
        for ( i=0; i<surface.nvoxels; i++)
        {
            for ( a=0; a<nactive; a++ ) {
                distance = get_distance ( &(surface.voxels[i].coords[0]), &(atoms[actatoms[a]].coords[0]) );
                surface.voxels[i].data += prefactor * exp( sqr( distance ) / (mttsqzeta));

            }
        }
    }
#endif

#ifndef OPTSURF
/* FU| check here, this is the plain code to do it
 */

    if ( periodic ) {
#ifdef OPENMP
#pragma omp parallel for default(none) \
    private(i,j,a,distance,skip) shared(surface,atoms,natoms,mask,prefactor,mttsqzeta,pbc,maxdist,cubedge,periodic) \
    // schedule(dynamic)
#endif
        for ( i=0; i<surface.nvoxels; i++)
        {
            for ( a=0; a<natoms; a++ ) {
                distance = get_distance_periodic ( &(surface.voxels[i].coords[0]), &(atoms[mask[a]].coords[0]), pbc );
                surface.voxels[i].data += prefactor * exp( sqr( distance ) / (mttsqzeta));
                }
        }
    }
    else {
#ifdef OPENMP
#pragma omp parallel for default(none) \
    private(i,j,a,distance,skip) shared(surface,atoms,nactive,actatoms,prefactor,mttsqzeta,pbc,maxdist,cubedge,periodic) \
    // schedule(dynamic)
#endif
        for ( i=0; i<surface.nvoxels; i++)
        {
            for ( a=0; a<nactive; a++ ) {
                distance = get_distance ( &(surface.voxels[i].coords[0]), &(atoms[actatoms[a]].coords[0]) );
                surface.voxels[i].data += prefactor * exp( sqr( distance ) / (mttsqzeta));

            }
        }
    }
#endif

    if ( !( periodic ) )
        free ( actatoms );

    char buf[MAXSTRLEN];

    if ( output )
    {
        sprintf(buf, "%s%s", outputprefix, "instant-surface.cube");
        write_cubefile(buf, &surface);
    }

    bindices = cubes_larger_than ( surfcut, &surface);
    sindices = invert_indices ( surface.nvoxels, bindices );

    if ( provide_mask ) {
        i = 0;
        while ( bindices[i] != -1 ) {
            surface.voxels[bindices[i]].data = ZERO;
            i++;
        }
    }

    i = 0;
    int nsurf = 0;
    while ( sindices[i] != -1 ) {
        if ( provide_mask )
            surface.voxels[sindices[i]].data = ONE;
        i++;
        nsurf = i;
    }

#if DEBUG
    printf("%i voxels belong to surface and occupy a volume of %21.10f Bohr^3\n", nsurf, nsurf*surface.dv);
#endif

    if ( output > 1)
    {
        sprintf(buf, "%s%s", outputprefix, "instant-surface-plain.cube");
        write_cubefile(buf, &surface);
    }

    free ( sindices );
    free ( bindices );
    return surface;
}

void get_2d_representation_ils ( cube_t * surface, real ** surf_2d_up, real ** surf_2d_down, real surfcut, int newsurf, int * surf_up_inds, int * surf_down_inds, int direction )
{
    int i, j, k;
    int baseind;
    int upper, lower;
    int cnt_up = 0;
    int cnt_down = 0;
    /* check here, and pass these values in function call or so */
    /* for now take a tenth of the whole surface-direction */
    int nproint = surface->n[direction] / 10;

    /* then interpolate to a hundresth of that spacing */
    int ninterp = nproint * 100;

    real * proint;
    real * interp;
    real * xvals;
    real * dx;

    dx = get_box_volels(surface);
#ifdef DEBUG
    printf("volume element: %f\n", dx[direction]);
#endif

    // real cutoff = surfcut * .9;
    /* check here, that we can do that, everything that is not surface should be zero'd out */
    /* check here, why does 0.5 work?, for interpolation, should we not use something that is not zero'd out? */

    xvals = ( real * ) malloc ( nproint * sizeof ( real ) );
    proint = ( real * ) malloc ( nproint * sizeof ( real ) );
    interp = ( real * ) malloc ( ninterp * sizeof ( real ) );

    if ( newsurf ) {
        surf_up_inds = ( int * ) malloc ( ( surface->nvoxels + 1 ) * sizeof ( int ) );
        surf_down_inds = ( int * ) malloc ( ( surface->nvoxels + 1 ) * sizeof ( int ) );
    }

// REINCLUDE
    // printf("%i %i\n", nproint, ninterp );

    /*
     * real ** surf_2d_up;
     * real ** surf_2d_down;
     */

    int crrnt, fndsrf;
    real tmpdt[DIM];
    real t, tfin;

    int surfup, surfdown;

    for ( i=0; i<surface->n[0]; i++ )
        for ( j=0; j<surface->n[1]; j++ ) {
            surfup = 0;
            surfdown = 0;
            for ( k=0; k<surface->n[2]; k++ ) {

                fndsrf = 0;

                baseind = surface->n[2] * ( j + surface->n[1] * i);

                // this is just the current voxel

                crrnt = baseind + k;

                /* check here, we still need to worry about the correct pbc treatment and exclude very first/last voxels */
                /* just simply wrap around boundaries, otherwise baseind -1 and +1 should provide correct boundary treament!? */


                // linear interpolation to find surface point:

                /* check here, new code should be done this way: */
                /* remember three voxels, one up and one down
                 * zero them once above the triple-loop
                 * zero them once a surface part was found
                 * if periodic make sure we treat boundaries correctly
                 * check if we find that we cross our pre-set cutoff (one (or two) voxel(s) > 0 and two (or one) voxels smaller than zero (or one equal to zero, only need to check adjacent points)
                 * if so interpolate between the three voxels and get surface point 
                 * save surface point in 1D array of indices
                 * if last voxel was found then add a marker for the end
                 * if more than two voxels were found then save the two extrema as surface points and the rest in a comment
                 * dunno right now how to handle the latter point exactly, because we would need to save it somehow in the loop, but we could
                 *     use a steadily growing 3D array that would carry the extra points
                 * and this, my friends is how we're going to do it
                 */

                /* check if we maybe get multiple or no assignments for either of the representations */
                /* check here, if we want the surface to be in arbitrary direction */
                /* do interpolation before assignment; it is the same for upper and lower surface */

                upper = crrnt + 1;
                lower = crrnt - 1;

                if ( k == 0 )
                    lower = baseind + surface->n[2] - 1;
                else if ( k == surface->n[2] - 1 )
                    lower = baseind;

                // this is not right!!!

                // if ( baseind == 0 )
                //     lower = lstvx;
                // else if ( baseind == lstvx )
                //     upper = 0;

                tmpdt[0] = surface->voxels[lower].data - surfcut;
                tmpdt[1] = surface->voxels[crrnt].data - surfcut;
                tmpdt[2] = surface->voxels[upper].data - surfcut;
// #ifdef DEBUG
//                 printf("%i %i %i\n", lower, crrnt, upper);
//                 printf("baseind: %i, values: %20.14f %20.14f %20.14f\n", i*j, tmpdt[0], tmpdt[1], tmpdt[2]);
// #endif

                // there are two obvious ways of doing it here
                //      I) always reassign surfup and surfdown, overhang ignored
                //      II) or do not reassign surfdown and always reassign surfdown (overhang will be counted as surface)
                //
                //      currently version II is done
                //      if we want to change that, also move the initialization of surfup and surfdown into the z-loop (that did not work so far) --> check here
                //      and also the reset of surfup to 0 in the conditional statement for fndsrf below --> check here

                // also not sure if the code works for closed interior surfaces (around a solute)
                // definitely not for interior surface in a cluster system

                // is the statement below sufficient?
                // maybe it's not, I've seen artifacts for 215 water molecule systems
                // and check here, if the code works fine (it does not for our crummy example, I know that much

                // check here, exchanged surfup and surfdown below, must've done it the wrong way 'round before, but please double-check
                if ( ( tmpdt[0] > 0 ) && ( tmpdt[2] < 0 ) ) {
                    fndsrf = crrnt;
                    // surfup = 1;
                    surfdown++;

                    if ( tmpdt[1] >= 0 ) {
                        t = lerp_to_t ( tmpdt[1], tmpdt[2], ZERO );
                        tfin = t;
                    }

                    if ( tmpdt[1] < 0 ) {
                        t = lerp_to_t ( tmpdt[0], tmpdt[1], ZERO );
                        tfin = t - 1;
                    }
                }

                else if ( ( tmpdt[0] < 0 ) && ( tmpdt[2] > 0 ) ) {
                    fndsrf = crrnt;
                    // surfdown++;
                    surfup = 1;

                    if ( tmpdt[1] >= 0 ) {
                        t = lerp_to_t ( tmpdt[0], tmpdt[1], ZERO );
                        tfin = t - 1;
                    }

                    if ( tmpdt[1] < 0 ) {
                        t = lerp_to_t ( tmpdt[1], tmpdt[2], ZERO );
                        tfin = t;
                    }

                }

                if ( fndsrf ) {
// #ifdef DEBUG
//                     printf("t is %f AND tfin is %f\n", t, tfin);
// #endif
                    if ( surfdown == 1 ) {

                        surf_2d_down[i][j] = surface->voxels[crrnt].coords[direction] + tfin * dx[direction];

                        /* check here and insert vertical interpolation */
                        if ( newsurf ) {
                            surf_up_inds[cnt_up] = lower;
                            cnt_up++;
                        }

                    }
                    else if ( surfup )

                        surf_2d_up[i][j] = surface->voxels[crrnt].coords[direction] + tfin * dx[direction];

                        /* check here and insert vertical interpolation */
                        if ( newsurf ) {
                            surf_down_inds[cnt_down] = upper;
                            cnt_down++;
                        }

                        // we could eually put this statement at the beginning of the loop
                        surfup = 0;
                }

            }
        }

    if ( newsurf ) {
        surf_up_inds[cnt_up] = -1;
        surf_down_inds[cnt_down] = -1;
    }

    free ( xvals );
    free ( proint );
    free ( interp );
}

void get_distance_to_surface ( real * disthi, real * distlo, int * inthi, int * intlo, cube_t * surface, real ** surf_2d_up, real ** surf_2d_down, atom_t * atoms, int * refmask, int nref, int natoms, real * pbc, int output, char * opref, int direction, real surfcut )
{
    int k, l;
    real ** lodi;
    real ** updi;

    lodi = allocate_matrix_real_2d ( surface->n[0], surface->n[1] );
    updi = allocate_matrix_real_2d ( surface->n[0], surface->n[1] );

    real spos[DIM];
    real com[DIM];

    get_center_of_mass ( com, atoms, refmask, nref);

    /* check here, we can do it like that because the surface was initialized to be cubic with (a,0,0), (0,a,0) and (0,0,a) as box vectors */
    for ( k=0; k<surface->n[0]; k++ )
        for ( l=0; l<surface->n[1]; l++ ) {
            /*check here, if these are really the centers of the voxels...*/
            spos[0] = k * surface->boxv[0][0];
            spos[1] = l * surface->boxv[1][1];
            /* check here, if we have the same X, Y description */

            spos[2] = surf_2d_down[k][l];
            lodi[k][l] = get_distance_periodic ( spos, com, pbc );

            spos[2] = surf_2d_up[k][l];
            updi[k][l] = get_distance_periodic ( spos, com, pbc );
        }

    int minx, miny, minz;

#ifdef DEBUG
    FILE *fsxyzlo, *fsxyzup;

    char tmp[MAXSTRLEN];
    sprintf(tmp, "%s%s", opref, "surflo.xyz");
    fsxyzlo = fopen(&tmp[0], "a");;

    sprintf(tmp, "%s%s", opref, "surfup.xyz");
    fsxyzup = fopen(&tmp[0], "a");;
#endif

    minx = 0;
    miny = 0;

    // printf ( "%i %i %i %i\n", minx, miny, surface->n[0], surface->n[1]);
    *distlo = lodi[minx][miny];

    /* check here, how to generalize to different orientation */

    minz = surf_2d_down[minx][miny] / surface->boxv[2][2];
    // printf("TEST 1: %f\n", surf_2d_down[minx][miny]);
    *intlo = get_index ( surface->n, minx, miny, minz );
    // printf("%i %i %i\n", minx, miny, minz);
    // printf("TEST 2: %f\n", surface->voxels[*intlo].coords[2]);

#ifdef DEBUG
    fprintf ( fsxyzlo, "%i\n\n", 2+natoms );
    fprintf ( fsxyzlo, "S %21.10f %21.10f %21.10f\n", BOHR * minx * surface->boxv[0][0], BOHR * miny * surface->boxv[1][1], BOHR * surf_2d_down[minx][miny] );
    fprintf ( fsxyzlo, "P %21.10f %21.10f %21.10f\n", BOHR * com[0], BOHR * com[1], BOHR * com[2] );

    int m;

    for ( m=0; m<natoms; m++ )
        fprintf ( fsxyzlo, "%s %21.10f %21.10f %21.10f\n", atoms[m].symbol, BOHR * atoms[m].coords[0], BOHR * atoms[m].coords[1], BOHR * atoms[m].coords[2] );

    fclose ( fsxyzlo );
#endif

    find_minimum_2d_real (&minx, &miny, updi, k, l );
    *disthi = updi[minx][miny];

    minz = surf_2d_up[minx][miny] / surface->boxv[2][2];
    *inthi = get_index ( surface->n, minx, miny, minz );
    // printf("%i %i %i\n", minx, miny, minz);

#ifdef DEBUG
    fprintf ( fsxyzup, "%i\n\n", 2+natoms );
    fprintf ( fsxyzup, "S %21.10f %21.10f %21.10f\n", BOHR * minx * surface->boxv[0][0], BOHR * miny * surface->boxv[1][1], BOHR * surf_2d_up[minx][miny] );
    fprintf ( fsxyzup, "P %21.10f %21.10f %21.10f\n", BOHR * com[0], BOHR * com[1], BOHR * com[2] );

    for ( m=0; m<natoms; m++ )
        fprintf ( fsxyzup, "%s %21.10f %21.10f %21.10f\n", atoms[m].symbol, BOHR * atoms[m].coords[0], BOHR * atoms[m].coords[1], BOHR * atoms[m].coords[2] );

    fclose ( fsxyzup );
#endif
    // printf("Minimum distances are for upper: %21.10f and lower: %21.10f\n", *distlo, *disthi);

    free_matrix_real_2d ( lodi, surface->n[0] );
    free_matrix_real_2d ( updi, surface->n[0] );

    if ( output > 1 ) {
        char funame[MAXSTRLEN];
        sprintf( funame, "%s%s", opref, "2dsurf-up.dat");

        write_matrix_real_2d_to_file_w_cont_spacing ( funame, surf_2d_up, surface->n[0], surface->n[1], &(surface->origin[0]), surface->boxv[0][0], surface->boxv[1][1] );

        char fdname[MAXSTRLEN];
        sprintf( fdname, "%s%s", opref, "2dsurf-down.dat");

        write_matrix_real_2d_to_file_w_cont_spacing ( fdname, surf_2d_down, surface->n[0], surface->n[1], &(surface->origin[0]), surface->boxv[0][0], surface->boxv[1][1] );
    }
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
            // check here, if that is the correct striding
            // NO, IT IS NOT, WHAT THE FUCK ARE YOU DOING, DUDE?
            dx[i] += cube->boxv[i][j];
            // NOW THIS SHOULD BE RIGHT
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
            // check here, if that is the correct striding
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
            // check here, if that is the correct striding
            dx[i] += cube->boxv[i][j];

    return dx;

}

