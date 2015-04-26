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
#include "surf.h"
#include "molmanipul.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int tstart;
int tstop;

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


    /* we will create an orthogonal box according to periodic boundary conditions and resolution input */
    /* i.e., it works for now only with orthogonal cells */
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

    surface = initialize_cube(orig, boxv, n, atoms, inpnatoms);

    sqzeta = sqr(zeta);
    real mttsqzeta = -2. * sqzeta;
    dummy = 2. * PI * sqzeta;
    prefactor = 1. / dummy / (sqrt(dummy));

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
    }

    /* check here, this will only work if our pbc always starts at (0,0,0) */
    for ( i=0; i<DIM; i++ ) {
        // round down to not overestimate things here
        /* check here, not sure which one is the correct one */

        addspace[i] = lroundf ( ( pbc[i] - cubpbc[i] + surface.origin[i] ) / dx[i] );
    }

    real distance;

    int skip;

    int * actatoms;
    int nactive = 0;

    int cnt = 0;

    if ( !(periodic ) ) {
        actatoms = ( int * ) malloc ( ( natoms + 1 ) * sizeof ( int ) );

        for ( a=0; a<natoms; a++ ) {
            skip = 0;

            for ( j=0; j<DIM; j++ ) {
                if ( ( atoms[mask[a]].coords[j] < cubedge[j] ) || ( atoms[mask[a]].coords[j] > cubedge[DIM+j] ) ) {
                    skip = 1;
                    break;
                }
            }

            if ( !( skip ) ) {
                actatoms[cnt] = mask[a];
                cnt++;
            }
        }
        actatoms[cnt] = -1;
        nactive = cnt;
    }

#if DEBUG
    printf("%i active atoms for surface generation\n", cnt);
#endif

#ifndef OLDSURF

    // printf("We are using the optimized, but not debugged version of the code!!!\n", trplzt);

    real mxdst;
    int mxvox;
    int mx[DIM];
    int mn[DIM];
    int k;
    int wrki, wrkj, wrkk;
    int index[DIM];
    int tmpndx;
    real resarr[DIM];
    real tmpdst;
    real cutshft = prefactor * exp( sqr( trplzt ) / (mttsqzeta));

    for ( k=0; k<DIM; k++ )
        resarr[k] = resolution;

    mxdst = trplzt + resolution;
    mxvox = mxdst / resolution;

#ifdef OPTSURF
    if ( periodic ) {

        for ( a=0; a<natoms; a++ ) {

            // get index of voxel where atom is sitting
            //
            get_index_triple ( index, &(atoms[mask[a]].coords[0]), pbc, resarr, periodic );

            // then determine mni, mxi, ... and so on

            for ( k=0; k<DIM; k++ ) {
                mn[k] = index[k] - mxvox;
                mx[k] = index[k] + mxvox;
            }

            // printf("atom #%i\n", a);

            for ( i=mn[0]; i<mx[0]; i++ )
                for ( j=mn[1]; j<mx[1]; j++ )
                    for ( k=mn[2]; k<mx[2]; k++ ) {

                        periodify_indices ( &wrki, &(surface.n[0]), &i, 1 );
                        periodify_indices ( &wrkj, &(surface.n[1]), &j, 1 );
                        periodify_indices ( &wrkk, &(surface.n[2]), &k, 1 );

                        tmpndx = get_index ( surface.n, wrki, wrkj, wrkk );

                        // printf("%i %i %i %i\n", wrki, wrkj, wrkk, tmpndx);

                        distance = get_distance_periodic ( &(surface.voxels[tmpndx].coords[0]), &(atoms[mask[a]].coords[0]), pbc );

                        if ( distance > trplzt )
                            continue;

                        surface.voxels[tmpndx].data += prefactor * exp( sqr( distance ) / (mttsqzeta));
                        surface.voxels[tmpndx].data -= cutshft;

                    }
        }
    }

#else

    if ( periodic ) {
#ifdef OPENMP
#pragma omp parallel for default(none) \
    private(i,j,k,a,distance,skip,tmpdst) shared(surface,atoms,natoms,mask,prefactor,mttsqzeta,pbc,maxdist,cubedge,periodic,trplzt,cutshft) \
        schedule(guided, surface.n[2])
    // schedule(dynamic)
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

                surface.voxels[i].data += prefactor * exp( sqr( distance ) / (mttsqzeta));
                surface.voxels[i].data -= cutshft;

                }
        }
    }
    else {
#ifdef OPENMP
#pragma omp parallel for default(none) \
    private(i,j,a,distance,skip) shared(surface,atoms,nactive,actatoms,prefactor,mttsqzeta,pbc,maxdist,cubedge,periodic,cutshft,trplzt)
#endif
        for ( i=0; i<surface.nvoxels; i++)
        {
            for ( a=0; a<nactive; a++ ) {
                distance = get_distance ( &(surface.voxels[i].coords[0]), &(atoms[actatoms[a]].coords[0]) );

                if ( distance > trplzt )
                    continue;

                surface.voxels[i].data += prefactor * exp( sqr( distance ) / (mttsqzeta));
                surface.voxels[i].data -= cutshft;

            }
        }
    }
#endif
#endif

#ifdef OLDSURF
    if ( periodic ) {
#ifdef OPENMP
#pragma omp parallel for default(none) \
    private(i,j,a,distance,skip) shared(surface,atoms,natoms,mask,prefactor,mttsqzeta,pbc,maxdist,cubedge,periodic)
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
    private(i,j,a,distance,skip) shared(surface,atoms,nactive,actatoms,prefactor,mttsqzeta,pbc,maxdist,cubedge,periodic)
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

real ** get_2d_representation_ils ( int * nsurf, int ** drctn, cube_t * surface, real surfcut, int newsurf, int * surf_inds )
{
    int i, j, k, l;
    int baseind;
    int upper, lower;
    int cnt_up = 0;
    int cnt_down = 0;
    int direction;
    /* check here, and pass these values in function call or so */
    /* for now take a tenth of the whole surface-direction */

    /* then interpolate to a hundresth of that spacing */

    real * dx;
    real tmpsrf;

    dx = get_box_volels(surface);
#ifdef DEBUG
    printf("volume element: %f\n", dx[direction]);
#endif

    if ( newsurf )
        surf_inds = ( int * ) malloc ( ( surface->nvoxels + 1 ) * sizeof ( int ) );

    int crrnt, fndsrf;
    real tmpdt[DIM];
    real t, tfin;
    real **surfpts;
    real dblsrfct = 2. * surfcut;

    // we'll assume every point could be a point at the surface and sort out things later on
    // we'll make no division between upper and lower surface, because it's bullshit anyways

    *nsurf = 0;

    surfpts = (real **) malloc ( DIM * sizeof ( real *) );
    int *tmpdir = ( int * ) malloc ( surface->nvoxels * sizeof ( int ) );

    for ( i=0; i<DIM; i++ )
        surfpts[i] = (real *) malloc ( surface->nvoxels * sizeof ( real ) );

    int surfup, surfdown;

    for ( i=0; i<surface->n[0]; i++ )
        for ( j=0; j<surface->n[1]; j++ ) {

            baseind = surface->n[2] * ( j + surface->n[1] * i);
            for ( k=0; k<surface->n[2]; k++ ) {

                fndsrf = 0;

                // defines indices for z-search (this is fine)
                crrnt = baseind + k;

                upper = crrnt + 1;
                lower = crrnt - 1;

                if ( k == 0 )
                    lower = baseind + surface->n[2] - 1;
                else if ( k == ( surface->n[2] - 1 ) )
                    upper = baseind;

                tmpdt[0] = surface->voxels[lower].data;
                tmpdt[1] = surface->voxels[crrnt].data;
                tmpdt[2] = surface->voxels[upper].data;

                // this will get all the fricking surface points
                // if ( ( tmpdt[2] + tmpdt[0] ) < dblsrfct ) {
                //     fndsrf = 1;
                //     direction = 0;
                // }

                if ( ( ( tmpdt[2] > surfcut ) && ( tmpdt[1] < surfcut ) && ( tmpdt[0] < surfcut ) ) ||
                ( ( tmpdt[2] < surfcut ) && (tmpdt[1] < surfcut) &&  ( tmpdt[0] > surfcut ) ) ) {
                    fndsrf = 1;
                    direction = 2;
                }

                // defines indices for y-search
                if ( ! ( fndsrf ) ) {
                    crrnt = baseind + k;

                    upper = get_index ( surface->n, i, (j+1), k );
                    lower = get_index ( surface->n, i, (j-1), k );

                    if ( j == 0 )
                        lower = get_index ( surface->n, i, surface->n[1]-1, k );
                    else if ( j == ( surface->n[1] - 1 ) )
                        upper = get_index ( surface->n, i, 0, k );

                    tmpdt[0] = surface->voxels[lower].data;
                    tmpdt[1] = surface->voxels[crrnt].data;
                    tmpdt[2] = surface->voxels[upper].data;

                    // if ( ( tmpdt[2] + tmpdt[0] ) < dblsrfct ) {
                    //     fndsrf = 1;
                    //     direction = 1;
                    // }

                    if ( ( ( tmpdt[2] > surfcut ) && ( tmpdt[1] < surfcut ) && ( tmpdt[0] < surfcut ) ) ||
                    ( ( tmpdt[2] < surfcut ) && (tmpdt[1] < surfcut) &&  ( tmpdt[0] > surfcut ) ) ) {
                        fndsrf = 1;
                        direction = 1;
                    }
                }

                // defines indices for x-search
                if ( ! ( fndsrf ) ) {
                    crrnt = baseind + k;

                    upper = get_index ( surface->n, (i+1), j, k );
                    lower = get_index ( surface->n, (i-1), j, k );

                    if ( i == 0 )
                        lower = get_index ( surface->n, surface->n[0]-1, j, k );
                    else if ( i == ( surface->n[0] - 1 ) )
                        upper = get_index ( surface->n, 0, j, k );

                    tmpdt[0] = surface->voxels[lower].data;
                    tmpdt[1] = surface->voxels[crrnt].data;
                    tmpdt[2] = surface->voxels[upper].data;

                    // if ( ( tmpdt[2] + tmpdt[0] ) < dblsrfct ) {
                    //     fndsrf = 1;
                    //     direction = 2;
                    // }

                    if ( ( ( tmpdt[2] > surfcut ) && ( tmpdt[1] < surfcut ) && ( tmpdt[0] < surfcut ) ) ||
                    ( ( tmpdt[2] < surfcut ) && (tmpdt[1] < surfcut) &&  ( tmpdt[0] > surfcut ) ) ) {
                        fndsrf = 1;
                        direction = 0;
                    }
                }

                // carefully check here again the assignment of the surface direction
                // maybe use a tri-linear interpolation (but i don't think it's actually needed
                // and substitute the stupid zeros by surfcut, that should work equally well

                if ( fndsrf ) {

                    if ( ( tmpdt[0] > surfcut ) && ( tmpdt[2] < surfcut ) ) {
                        fndsrf = crrnt;
                        // surfup = 1;
                        surfup++;

                        // change the operators below (< and >=) check here for validity
                        if ( tmpdt[2] < surfcut ) {
                            t = lerp_to_t ( tmpdt[1], tmpdt[2], surfcut );
                            tfin = t;
                        }

                        if ( tmpdt[0] >= surfcut ) {
                            t = lerp_to_t ( tmpdt[0], tmpdt[1], surfcut );
                            tfin = t - 1;
                        }
                    }

                    else if ( ( tmpdt[0] < surfcut ) && ( tmpdt[2] > surfcut ) ) {
                        fndsrf = crrnt;
                        // surfdown++;
                        surfdown = 1;

                        if ( tmpdt[2] >= surfcut ) {
                            t = lerp_to_t ( tmpdt[0], tmpdt[1], surfcut );
                            tfin = t - 1;
                        }

                        if ( tmpdt[0] < surfcut ) {
                            t = lerp_to_t ( tmpdt[1], tmpdt[2], surfcut );
                            tfin = t;
                        }

                    }
                    // for now interpolate only in the direction that the surface is actually present
                    // printf("direction: %i\n", direction);
                    tmpsrf = surface->voxels[crrnt].coords[direction] + tfin * dx[direction];
                    tmpdir[*nsurf] = direction;

                    for ( l=0; l<DIM; l++ ) {
                        if ( l == direction )
                            surfpts[l][*nsurf] = tmpsrf;
                        else
                            surfpts[l][*nsurf] = surface->voxels[crrnt].coords[l];
                    }

                    (*nsurf)++;

                    // /* check here and insert vertical interpolation */
                    // if ( newsurf )
                    //     surf_inds[*nsurf] = lower;

                }

            }
        }

    // if ( newsurf )
    //     surf_inds[*nsurf] = -1;

    real ** finalsurf;

    // check here and insert allocate_matrix_real_2d (also above);
    finalsurf = ( real ** ) malloc ( *nsurf * sizeof ( real * ) );
    *drctn =  ( int * ) malloc ( *nsurf * sizeof ( int ) );

    for ( i=0; i<*nsurf; i++ ) {
        finalsurf[i] = ( real * ) malloc ( DIM * sizeof ( real ) );

        for ( k=0; k<DIM; k++ )
            finalsurf[i][k] = surfpts[k][i];
    }

    for ( i=0; i<DIM; i++ )
        free ( surfpts[i] );

    free ( surfpts );
    free ( tmpdir );

    free ( dx );

    return finalsurf;
}

real get_distance_to_surface ( cube_t * surface, int nsurf, real ** surfpts, int * direction, atom_t * atoms, int * refmask, int nref, int natoms, real * pbc, int output, char * opref, real surfcut )
{
    int k, l;
    real ** lodi;
    real ** updi;
    real * dx;

    dx = get_box_volels(surface);
    real * dsts = ( real * ) malloc ( nsurf * sizeof ( real ) );

    real spos[DIM];
    real com[DIM];
    real dstnc;

    get_center_of_mass ( com, atoms, refmask, nref);

    for ( k=0; k<nsurf; k++ )
        dsts[k] = get_distance_periodic ( surfpts[k], com, pbc );

#ifdef DEBUG
    FILE *fsxyzlo, *fsxyzup;

    char tmp[MAXSTRLEN];
    sprintf(tmp, "%s%s", opref, "surflo.xyz");
    fsxyzlo = fopen(&tmp[0], "a");;

    sprintf(tmp, "%s%s", opref, "surfup.xyz");
    fsxyzup = fopen(&tmp[0], "a");;
#endif

    int min = 0;

    dstnc = find_minimum_1d_real (&min, dsts, nsurf );

    // now it's more difficult to check whether in- or outside of bulk (because we don't know the direction of the surface at that point, unless we safe it maybe)
    // that's why we have the direction array now

    // if ( surf_2d_down[minx][miny] > com[2] )
    //     *distlo *= -1;

    free ( dsts );
    free ( dx );

    return dstnc;

//     if ( output > 1 ) {
//         char funame[MAXSTRLEN];
//         sprintf( funame, "%s%s", opref, "2dsurf-up.dat");
//
//         write_matrix_real_2d_to_file_w_cont_spacing ( funame, surf_2d_up, surface->n[0], surface->n[1], &(surface->origin[0]), surface->boxv[0][0], surface->boxv[1][1] );
//
//         char fdname[MAXSTRLEN];
//         sprintf( fdname, "%s%s", opref, "2dsurf-down.dat");
//
//         write_matrix_real_2d_to_file_w_cont_spacing ( fdname, surf_2d_down, surface->n[0], surface->n[1], &(surface->origin[0]), surface->boxv[0][0], surface->boxv[1][1] );
//     }
}
