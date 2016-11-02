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
#include "errors.h"
#include "molmanipul.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int tstart;
int tstop;

cube_t instant_surface_periodic ( int * mask, atom_t * atoms, int inpnatoms, real *zeta, real surfcut, int output, char * outputprefix, real * pbc, real resolution, real accuracy, int provide_box, real * origincube, int * ncube, real boxvcube[DIM][DIM], int periodic, int provide_mask )
{
    int i, j, natoms;
    int a;
    int * sindices;
    int * bindices;
    real sqzeta;
    cube_t surface;
    // real prefactor, dummy, cutshft;

    real orig[DIM];
    real boxv[DIM][DIM];
    real refc[DIM];
    int n[DIM];

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

    int wrki, wrkj, wrkk;
    int index[DIM];
    int tmpndx;
    real resarr[DIM];
    real tmpdst;

    int k;
    for ( k=0; k<DIM; k++ )
        resarr[k] = resolution;

#ifdef OPENMP
#pragma omp parallel for default(none) \
    private(a,i,j,k,wrki,wrkj,wrkk,index,tmpndx,distance) shared(atoms,pbc,resarr,periodic,surface,natoms,mask,zeta,resolution) // \
        // schedule(guided, surface.n[2])
    // schedule(dynamic)
#endif
    //this natoms here is already the one accounting for number of atoms in mask only
    for ( a=0; a<natoms; a++ ) {

        real sqzeta = sqr(zeta[mask[a]]);
        real trplzt = 3*zeta[mask[a]];

        real mttsqzeta = -2. * sqzeta;
        real dummy = 2. * PI * sqzeta;

        real prefactor = 1. / dummy / (sqrt(dummy));
        real cutshft = prefactor * exp( sqr( trplzt ) / (mttsqzeta));

        real mxdst = trplzt + 2*resolution;
        real mxvox[DIM];

        int mx[DIM];
        int mn[DIM];
        int k;

        for ( k=0 ; k<DIM; k++ ) {
            mxvox[k] = mxdst / resolution;

            if ( 2*mxvox[k] > surface.n[k] )
                mxvox[k] -= (int) ceil ( (real) (2*mxvox[k] - surface.n[k]) / 2. );
        }

        // this should not happen anymore
        for ( k=0; k<DIM; k++ ) {
            if ( ( 2*mxvox[k] ) > ( surface.n[k] ) ) {
                print_error ( PROGRAM_BROKEN, "to use box dimension smaller than two times surface calculation cutoff");
                exit ( PROGRAM_BROKEN );
            }
        }

        // get index of voxel where atom is sitting
        //
        get_index_triple ( index, &(atoms[mask[a]].coords[0]), pbc, surface.origin, surface.n, resarr, periodic );

        // then determine mni, mxi, ... and so on

        for ( k=0; k<DIM; k++ ) {
            mn[k] = index[k] - mxvox[k];
            mx[k] = index[k] + mxvox[k];
        }

        for ( i=mn[0]; i<mx[0]; i++ )
            for ( j=mn[1]; j<mx[1]; j++ )
                for ( k=mn[2]; k<mx[2]; k++ ) {

                    periodify_indices ( &wrki, &(surface.n[0]), &i, 1 );
                    periodify_indices ( &wrkj, &(surface.n[1]), &j, 1 );
                    periodify_indices ( &wrkk, &(surface.n[2]), &k, 1 );

                    tmpndx = get_index ( surface.n, wrki, wrkj, wrkk );

                    if ( periodic )
                        distance = get_distance_periodic ( &(surface.voxels[tmpndx].coords[0]), &(atoms[mask[a]].coords[0]), pbc );
                    else
                        distance = get_distance ( &(surface.voxels[tmpndx].coords[0]), &(atoms[mask[a]].coords[0]) );

                    if ( distance > trplzt )
                        continue;

                    // the below formula can be simplified, check here
#pragma omp atomic update
                    surface.voxels[tmpndx].data += prefactor * exp( sqr( distance ) / (mttsqzeta)) - cutshft;

                }
    }

    if ( !( periodic ) )
        free ( actatoms );

    char buf[MAXSTRLEN];

    if ( output > 1 )
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

    if ( ( provide_mask ) && ( output > 2) ) {
        sprintf(buf, "%s%s", outputprefix, "instant-surface-plain.cube");
        write_cubefile(buf, &surface);
    }

    free ( sindices );
    free ( bindices );
    return surface;
}

real ** get_2d_representation_ils ( int * nsurf, int ** drctn, real ** grad, cube_t * surface, real surfcut, int newsurf, int * surf_inds, int direction, real * area, int periodic )
{
    int i, j, k, l;
    int baseind;
    int upper, lower;
    int cnt_up = 0;
    int cnt_down = 0;

    real * dx;
    real tmpsrf;

    int d;
    int ix[DIM];

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

    int  srf[DIM];
    int lrdfnd;

    real grd[DIM];
    real tpx[DIM];
    real dlt;

    // we'll assume every point could be a point at the surface and sort out things later on
    // we'll make no division between upper and lower surface, because it's bullshit anyways

    *nsurf = 0;

    surfpts = (real **) malloc ( DIM * sizeof ( real *) );
    int *tmpdir = ( int * ) malloc ( surface->nvoxels * sizeof ( int ) );
    real *tgrd = ( real * ) malloc ( surface->nvoxels * sizeof ( real ) );

    real mxgrd;
    real dt[DIM];
    real tmparea[DIM];

    int lofin, hifin;
    int dstrt, dstp;

    if ( ( direction >= 0 ) && ( direction < DIM ) ) {
        dstrt = direction;
        dstp = direction + 1;
    }
    else {
        dstrt = 0;
        dstp = DIM;
    }

    get_box_areas_pointer (tmparea, surface, dx );

    for ( i=0; i<DIM; i++ )
        surfpts[i] = (real *) malloc ( surface->nvoxels * sizeof ( real ) );

    int mxd[DIM];
    int mnd[DIM];

    if ( periodic )
        for ( k=0; k<DIM; k++ ) {
            mxd[k] = surface->n[k];
            mnd[k] = 0;
        }
    else
        for ( k=0; k<DIM; k++ ) {
            mxd[k] = surface->n[k] - 1;
            mnd[k] = 1;
        }

    for ( i=mnd[0]; i<mxd[0]; i++ )
        for ( j=mnd[1]; j<mxd[1]; j++ ) {

            baseind = surface->n[2] * ( j + surface->n[1] * i);
            for ( k=mnd[2]; k<mxd[2]; k++ ) {

                crrnt = baseind + k;
                fndsrf = 0;
                lrdfnd = 0;

                ix[0] = i;
                ix[1] = j;
                ix[2] = k;

                mxgrd = ZERO;

                for ( d=dstrt; d<dstp; d++ ) {

                    fndsrf = check_if_surface_voxel ( &upper, &lower, tmpdt, surface, ix, d, surfcut, periodic );

                    tpx[d] = ZERO;
                    lrdfnd += fndsrf;
                    srf[d] = fndsrf;

                    if ( fndsrf ) {
                        grd[d] = surface->voxels[upper].data - surface->voxels[lower].data;
                        grd[d] /= 2 * dx[d];

                        if ( ( fabs ( grd[d] ) ) > ( fabs ( mxgrd ) ) ) {
                            tmpdir[*nsurf] = d;
                            mxgrd = grd[d];
                            tgrd[*nsurf] = mxgrd;

                            for ( l=0; l<DIM; l++ )
                                dt[l] = tmpdt[l];
                        }
                    }

                }

                if ( lrdfnd ) {

                    for ( d=0; d<DIM; d++ ) {
                        if ( d == tmpdir[*nsurf] ) {

                            *area += tmparea[d];

                            t = lerp_to_t ( dt[1], dt[2], surfcut );

                            tpx[d] = surface->voxels[crrnt].coords[d] + t * dx[d];

                        }
                        else {
                            tpx[d] = surface->voxels[crrnt].coords[d];
                        }

                        surfpts[d][*nsurf] = tpx[d];

                    }

                (*nsurf)++;

                }

            }
        }

    if ( newsurf )
        surf_inds[*nsurf] = -1;

    real ** finalsurf;

    // check here and insert allocate_matrix_real_2d (also above);
    finalsurf = ( real ** ) malloc ( *nsurf * sizeof ( real * ) );
    *drctn =  ( int * ) malloc ( *nsurf * sizeof ( int ) );
    *grad = ( real * ) malloc ( *nsurf * sizeof ( real ) );

    for ( i=0; i<*nsurf; i++ ) {

        (*drctn)[i] = tmpdir[i];
        (*grad)[i] = tgrd[i];

        finalsurf[i] = ( real * ) malloc ( DIM * sizeof ( real ) );
        for ( k=0; k<DIM; k++ )
            finalsurf[i][k] = surfpts[k][i];
    }

    for ( i=0; i<DIM; i++ )
        free ( surfpts[i] );

    free ( surfpts );
    free ( tmpdir );
    free ( tgrd );

    free ( dx );

    return finalsurf;
}

real get_distance_to_surface ( int * mnnd, cube_t * surface, int nsurf, real ** surfpts, int * direction, real * grad, atom_t * atoms, int * refmask, int nref, int natoms, real * pbc, int output, char * opref, real surfcut, int periodic )
{
    int k, l;
    int crrnt[DIM];
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
        if ( periodic )
            dsts[k] = get_distance_periodic ( surfpts[k], com, pbc );
        else
            dsts[k] = get_distance ( surfpts[k], com );

    int min = 0;

    dstnc = find_minimum_1d_real (&min, dsts, nsurf );

    get_index_triple ( crrnt, com, pbc, surface->origin, surface->n, dx, periodic );

    int lower[DIM], upper[DIM];

    for ( l=0; l<DIM; l++ ) {
        if ( l == direction[min] ) {
            lower[l] = crrnt[l] - 1;
            upper[l] = crrnt[l] + 1;
        }
        else {
            lower[l] = crrnt[l];
            upper[l] = crrnt[l];
        }
    }

    int loind, hiind;

    periodify_indices ( lower, surface->n, lower, DIM );
    periodify_indices ( upper, surface->n, upper, DIM );

    loind = get_index ( surface->n, lower[0], lower[1], lower[2] );
    hiind = get_index ( surface->n, upper[0], upper[1], upper[2] );

    if ( ( ( grad[min] < 0. ) && ( com[direction[min]] > surfpts[min][direction[min]] ) ) ||
         ( ( grad[min] > 0. ) && ( com[direction[min]] < surfpts[min][direction[min]] ) ) ) {
        // for later checks if solute is above/below surface
        // printf("%5i\n", refmask[0]);
        dstnc *= -1.;
    }

    free ( dsts );
    free ( dx );

    *mnnd = min;
    return dstnc;

}

int check_if_surface_voxel ( int * upper, int * lower, real * tmpdt, cube_t * surface, int * ix, int direction, real surfcut, int periodic )
{
    int k, l, ihi, ilo;
    signed int inds[DIM][DIM];
    int voxinds[DIM];
    int fndsrf = 0;

    for ( k=0; k<DIM; k++ )
        for ( l=0; l<DIM; l++ )
            inds[k][l] = ix[l];

    // check here, this will invert the gradient calculation
    inds[0][direction] -= 1;
    inds[2][direction] += 1;

    if ( periodic ) {
        periodify_indices ( &(inds[0][direction]), &(surface->n[direction]), &(inds[0][direction]), 1 );
        periodify_indices ( &(inds[2][direction]), &(surface->n[direction]), &(inds[2][direction]), 1 );
    }

    for ( k=0; k<DIM; k++ ) {
        voxinds[k] = get_index ( surface->n, inds[k][0], inds[k][1], inds[k][2] );
        tmpdt[k] = surface->voxels[voxinds[k]].data;
    }

#ifdef ALLSURFPOINTS
    if ( ( tmpdt[2] + tmpdt[0] ) < 2 * surfcut ) {
        fndsrf = 1;
    }
#else
    // or maybe the other way 'round
    // we could also set fndsrf to something more informative that will tell us which two voxels we should use!?
    // the way below, the middle voxels will always be close to the desired surface segment!?
    // if ( ( ( tmpdt[2] < surfcut ) && ( tmpdt[1] > surfcut ) && ( tmpdt[0] > surfcut ) ) ||
    // ( ( tmpdt[2] > surfcut ) && (tmpdt[1] > surfcut) &&  ( tmpdt[0] < surfcut ) ) )

    // why exactly does this not work?
    // if ( ( ( tmpdt[1] > surfcut ) && ( tmpdt[2] < surfcut ) ) ||
    //      ( ( tmpdt[1] < surfcut ) && ( tmpdt[2] > surfcut ) ) )

    if ( ( ( tmpdt[2] > surfcut ) && ( tmpdt[1] < surfcut ) && ( tmpdt[0] < surfcut ) ) ||
         ( ( tmpdt[2] < surfcut ) && ( tmpdt[1] < surfcut ) && ( tmpdt[0] > surfcut ) ) )
        fndsrf = 1;
#endif

    *upper = voxinds[2];
    *lower = voxinds[0];

    return fndsrf;
}

real get_bulk_volume ( cube_t * surface, real surfcut )
{
    int i;
    real vol = ZERO;

    for ( i=0; i<surface->nvoxels; i++ )
        if ( surface->voxels[i].data >= surfcut )
            vol += ONE;

    vol *= surface->dv;

    return vol;
}
