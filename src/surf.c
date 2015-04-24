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

                        if ( i < 0 )
                            wrki = i + surface.n[0];
                        else if ( i >= surface.n[0] )
                            wrki = i - surface.n[0];
                        else
                            wrki = i;

                        if ( j < 0 )
                            wrkj = j + surface.n[1];
                        else if ( j >= surface.n[1] )
                            wrkj = j - surface.n[1];
                        else
                            wrkj = j;

                        if ( k < 0 )
                            wrkk = k + surface.n[2];
                        else if ( k >= surface.n[2] )
                            wrkk = k - surface.n[2];
                        else
                            wrkk = k;

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

    xvals = ( real * ) malloc ( nproint * sizeof ( real ) );
    proint = ( real * ) malloc ( nproint * sizeof ( real ) );
    interp = ( real * ) malloc ( ninterp * sizeof ( real ) );

    if ( newsurf ) {
        surf_up_inds = ( int * ) malloc ( ( surface->nvoxels + 1 ) * sizeof ( int ) );
        surf_down_inds = ( int * ) malloc ( ( surface->nvoxels + 1 ) * sizeof ( int ) );
    }

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

                crrnt = baseind + k;

                upper = crrnt + 1;
                lower = crrnt - 1;

                if ( k == 0 )
                    lower = baseind + surface->n[2] - 1;
                else if ( k == surface->n[2] - 1 )
                    upper = baseind;

                tmpdt[0] = surface->voxels[lower].data;
                tmpdt[1] = surface->voxels[crrnt].data;
                tmpdt[2] = surface->voxels[upper].data;

                // carefully check here again the assignment of the surface direction
                // maybe use a tri-linear interpolation (but i don't think it's actually needed
                // and substitute the stupid zeros by surfcut, that should work equally well

                if ( ( tmpdt[0] > surfcut ) && ( tmpdt[2] < surfcut ) ) {
                    fndsrf = crrnt;
                    // surfup = 1;
                    surfdown++;

                    if ( tmpdt[1] >= surfcut ) {
                        t = lerp_to_t ( tmpdt[1], tmpdt[2], surfcut );
                        tfin = t;
                    }

                    if ( tmpdt[1] < surfcut ) {
                        t = lerp_to_t ( tmpdt[0], tmpdt[1], surfcut );
                        tfin = t - 1;
                    }
                }

                else if ( ( tmpdt[0] < surfcut ) && ( tmpdt[2] > surfcut ) ) {
                    fndsrf = crrnt;
                    // surfdown++;
                    surfup = 1;

                    if ( tmpdt[1] >= surfcut ) {
                        t = lerp_to_t ( tmpdt[0], tmpdt[1], surfcut );
                        tfin = t - 1;
                    }

                    if ( tmpdt[1] < surfcut ) {
                        t = lerp_to_t ( tmpdt[1], tmpdt[2], surfcut );
                        tfin = t;
                    }

                }

                if ( fndsrf ) {
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

    free ( dx );
    free ( xvals );
    free ( proint );
    free ( interp );
}

void get_distance_to_surface ( real * disthi, real * distlo, int * inthi, int * intlo, cube_t * surface, real ** surf_2d_up, real ** surf_2d_down, atom_t * atoms, int * refmask, int nref, int natoms, real * pbc, int output, char * opref, int direction, real surfcut )
{
    int k, l;
    real ** lodi;
    real ** updi;
    real * dx;

    dx = get_box_volels(surface);
    lodi = allocate_matrix_real_2d ( surface->n[0], surface->n[1] );
    updi = allocate_matrix_real_2d ( surface->n[0], surface->n[1] );

    real spos[DIM];
    real com[DIM];

    get_center_of_mass ( com, atoms, refmask, nref);

    for ( k=0; k<surface->n[0]; k++ )
        for ( l=0; l<surface->n[1]; l++ ) {
            /*check here, if these are really the centers of the voxels...*/
            /*yes, please do that*/

            spos[0] = k * dx[0];
            spos[1] = l * dx[1];

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

    find_minimum_2d_real (&minx, &miny, lodi, surface->n[0], surface->n[1] );
    *distlo = lodi[minx][miny];

    // maybe make the below dependent on whether we are above/below the surface
    minz = (int) floor ( ( surf_2d_down[minx][miny] - surface->origin[2] ) / dx[2] );
    *intlo = get_index ( surface->n, minx, miny, minz );

    if ( surf_2d_down[minx][miny] < com[2] )
        *distlo *= -1;

    find_minimum_2d_real (&minx, &miny, updi, surface->n[0], surface->n[1] );
    *disthi = updi[minx][miny];

    minz = (int) ceil ( ( surf_2d_up[minx][miny] - surface->origin[2] ) / dx[2] );
    *inthi = get_index ( surface->n, minx, miny, minz );

    // check here for possible mishabs that could happen (not sure, there might be some spurious cases here (check here))
    if ( surf_2d_up[minx][miny] > com[2] )
        *disthi *= -1;

#ifdef DEBUG
    fprintf ( fsxyzlo, "%i\n\n", 2+natoms );
    fprintf ( fsxyzlo, "S %21.10f %21.10f %21.10f\n", BOHR * minx * surface->boxv[0][0], BOHR * miny * surface->boxv[1][1], BOHR * surf_2d_down[minx][miny] );
    fprintf ( fsxyzlo, "P %21.10f %21.10f %21.10f\n", BOHR * com[0], BOHR * com[1], BOHR * com[2] );

    int m;

    for ( m=0; m<natoms; m++ )
        fprintf ( fsxyzlo, "%s %21.10f %21.10f %21.10f\n", atoms[m].symbol, BOHR * atoms[m].coords[0], BOHR * atoms[m].coords[1], BOHR * atoms[m].coords[2] );

    fclose ( fsxyzlo );

    fprintf ( fsxyzup, "%i\n\n", 2+natoms );
    fprintf ( fsxyzup, "S %21.10f %21.10f %21.10f\n", BOHR * minx * surface->boxv[0][0], BOHR * miny * surface->boxv[1][1], BOHR * surf_2d_up[minx][miny] );
    fprintf ( fsxyzup, "P %21.10f %21.10f %21.10f\n", BOHR * com[0], BOHR * com[1], BOHR * com[2] );

    for ( m=0; m<natoms; m++ )
        fprintf ( fsxyzup, "%s %21.10f %21.10f %21.10f\n", atoms[m].symbol, BOHR * atoms[m].coords[0], BOHR * atoms[m].coords[1], BOHR * atoms[m].coords[2] );

    fclose ( fsxyzup );
#endif

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

    free ( dx );
}
