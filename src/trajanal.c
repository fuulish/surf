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
#include "trajanal.h"
#include "time.h"
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "errors.h"

int tstart;
int tstop;

int tanalize ( input_t * inppar )
{
    int snapsize;
    XDRFILE * xd_read;
    FILE * fxmol;
    // FILE * frepxyz;
    FILE * fsdist;
    // FILE * fxyz;
#ifdef DEBUG
    // FILE * fsxyzlo;
    // FILE * fsxyzup;
#endif

    int i;
    // int n[DIM];
    int nref;
    int natoms;
    int nmask;

    int * mask;
    int * refmask;

    atom_t * atoms;

    i = 0;
    nref = 0;
    nmask = 0;

    if ( !( inppar->periodic ) )
        printf("W A R N I N G: calculation proceeds w/o periodic boundary conditions. Often this is not what you want and Frank always forgot to set 'periodic' in input\n");

    /* read initial snapshot to get structure */

    char text[MAXSTRLEN];
    
    if ( inppar->xdrread ) {
        int result_xtc;
        int natoms_xtc;

        fxmol = fopen(&inppar->structure[0], "r");
        natoms = atoi( fgets ( text, MAXSTRLEN, fxmol ) );
        read_xmol(&inppar->structure[0], &atoms);

#ifdef XDRCAP
        xd_read = xdrfile_open ( &inppar->trajectory[0], "r" );
        result_xtc = read_xtc_natoms( &inppar->trajectory[0], &natoms_xtc);
#endif

        if ( exdrOK != result_xtc ) {
            printf("Something went wrong opening XDR file\n"); // Error
            exit ( 1 );
        }

        if ( natoms_xtc != natoms ) {
            printf("XDR file and xyz file do not contain the same number of atoms\nNot continuing\n");
            exit ( 1 );
        }

        fclose ( fxmol );
            
    }
    else {
        fxmol = fopen(&inppar->trajectory[0], "r");
        natoms = atoi( fgets ( text, MAXSTRLEN, fxmol ) );
        snapsize = xmol_snap_bytesize(fxmol);
        read_xmol(inppar->trajectory, &atoms);
    }

    mask = get_mask(inppar->maskkind, inppar->mask, inppar->nkinds, atoms, natoms);

    if ( strstr ( inppar->refmask, EMPTY ) == NULL )
        refmask = get_mask(inppar->refmaskkind, inppar->refmask, inppar->refnkinds, atoms, natoms);
    else {
        printf("Cannot do analysis without 'refmask'\n");
        exit ( 1 );
    }

    nref = 0;
    while ( refmask[nref] != -1 )
        nref++;

    if ( !(nref) ) {
        printf("Cannot continue with 0 reference atoms\n");
        exit ( 1 );
    }

    atom_t refatoms[nref];

    nmask = 0;
    while ( mask[nmask] != -1 )
        nmask++;

    for ( i=0; i<nref; i++ )
        refatoms[i] = atoms[refmask[i]];

    real *densprof;
    real mxdim;
    int ndprof;
    real drdprof = inppar->profileres;

    if ( inppar->tasknum == SURFDENSPROF ) {

        int mini;

        if ( ! ( inppar->pbcset ) ) {
             print_error ( MISSING_INPUT_PARAM, "pbc" );
             // check here, need to introduce return value for tanalize
             return MISSING_INPUT_PARAM;
        }

        // mxdim = find_maximum_1d_real ( &mini, inppar->pbc, DIM );
        // just use the information given by the 'direction' keyword

        mxdim = ZERO;
        for ( i=0; i<DIM; i++ )
            mxdim += sqr ( inppar->pbc[i] );

        mxdim = sqrt ( mxdim );

        // we want to use both above and below surface
        ndprof = 2 * ( int ) ( mxdim / inppar->profileres );

        densprof = ( real * ) calloc ( ndprof, sizeof(real) );

    }

    // real origin[DIM];
    // real boxv[DIM][DIM];

    // for (i=0; i<DIM; i++) {
    //     origin[i] = ZERO;
    //     n[i] = inppar->pbc[i] / inppar->resolution;

    //     for ( j=0; j<DIM; j++ )
    //         boxv[i][j] = ZERO;

    //     boxv[i][i] = inppar->pbc[i] / n[i];
    // }

    int counter = 0;

    int ntot = ( inppar->stop - inppar->start ) / inppar->stride;
    int frwrd = inppar->start + 1;
    char * htw = "w";

#ifdef OPTSURF
    if ( ( inppar->tasknum == SURFDIST ) || ( inppar->tasknum == SURFDENSPROF ) )
        printf("Using the optimized surface routine.\nIt is about 20+ times faster, but not debugged yet.\nAll tests so far give identical numerical results to older, naive version!\n\n");
#endif

    for ( i=inppar->start; i<inppar->stop; i += inppar->stride )
    {
        /* read snapshot */
        if ( inppar->xdrread ) {
            read_xtr_forward ( xd_read, frwrd, atoms, natoms );
            frwrd = inppar->stride;
        }
        else {
            xmolreader(fxmol, snapsize, i, atoms, natoms);
        }

        if ( ( inppar->tasknum == SURFDIST ) || ( inppar->tasknum == SURFDENSPROF ) ) {

            // FU| check here, still need surface area determination

            cube_t surface;
            int fake_n[DIM];
            real fake_origin[DIM];
            real fake_boxv[DIM][DIM];
            real disthi, distlo;
            int * surf_inds;
            int inthi, intlo;
            char tmp[MAXSTRLEN];

            char opref[MAXSTRLEN];
            sprintf(opref, "%s%i_", inppar->outputprefix, i);

            surface = instant_surface_periodic ( mask, atoms, natoms, inppar->zeta, inppar->surfacecutoff, inppar->output, opref, inppar->pbc, inppar->resolution, inppar->accuracy, 0, fake_origin, fake_n, fake_boxv, inppar->periodic, 0 );

            /* check here, and remove hard-coded surface direction */
            int * direction;
            int newsurf = 0;
            int nsurf = 0;

            real ** surfpts;

            surfpts = get_2d_representation_ils ( &nsurf, &direction, &surface, inppar->surfacecutoff, newsurf, surf_inds );

            // use function write_combined_xmol
            if ( inppar->surfxyz ) {
                FILE *fsxyzal;

                sprintf(tmp, "%s%i_%s", inppar->outputprefix, i, "atrep_surface.xyz");
                fsxyzal = fopen(&tmp[0], "w");;

                int a, g, k;

                fprintf ( fsxyzal, "%i\n\n", surface.natoms+nsurf );

                for ( a=0; a<surface.natoms; a++ ) {
                    fprintf ( fsxyzal, "    %s", surface.atoms[a].symbol );
                    for ( k=0; k<DIM; k++ ) {
                        fprintf ( fsxyzal, "    %21.10f", surface.atoms[a].coords[k]*BOHR );
                    }

                    fprintf ( fsxyzal, "\n");
                }

                for ( g=0; g<nsurf; g++ ) {
                    fprintf(fsxyzal, "%5s", "S");
                    for ( k=0; k<DIM; k++ ) {
                        fprintf ( fsxyzal, "    %21.10f", BOHR * surfpts[g][k]);
                    }
                    fprintf( fsxyzal, "\n" );
                }

                fclose ( fsxyzal );
            }

            real dstnc;
            if ( inppar->tasknum == SURFDENSPROF ) {

                int r;
                signed int ind;
                int hndprof = ndprof / 2.;

                for ( r=0; r<nref; r++ ) {

                    dstnc = get_distance_to_surface ( &surface, nsurf, surfpts, direction, atoms, &(refmask[r]), 1, natoms, inppar->pbc, inppar->output, opref, inppar->surfacecutoff, inppar->periodic );

                    ind = ( int ) floor ( dstnc / inppar->profileres );

                    // periodify_indices ( &ind, &ndprof, &ind, 1 );

                    densprof[ hndprof + ind ] += 1.;

                }
            }
            else {
                // check here, this needs to be done for all of the solute atoms, not just assume that there is only one
                dstnc = get_distance_to_surface ( &surface, nsurf, surfpts, direction, atoms, refmask, nref, natoms, inppar->pbc, inppar->output, opref, inppar->surfacecutoff, inppar->periodic );
            }

            /* check here, and move stuff for refinement box creation somewhere else */
        
            sprintf(tmp, "%s%s", inppar->outputprefix, "surfdist.dat");
            fsdist = fopen(&tmp[0], htw);

            if ( strncmp ( htw, "w", 1 ) == 0 ) {
                fprintf ( fsdist, "#            lower                 upper\n");
            }

            fprintf ( fsdist, "%21.10f %21.10f\n", distlo, disthi);
            fclose ( fsdist );
            htw = "a";

            int k;
            for ( k=0; k<nsurf; k++ )
                free ( surfpts[k] );

            free ( surfpts );
            free ( direction );
            free ( surface.atoms );
            free ( surface.voxels );
        }

        printf("%4.2f %% done\r", (real) counter / ntot * 100.);
        fflush(stdout);

        counter++;

    }
    printf("%4.2f %% done\n", (real) counter / ntot * 100.);

    if ( ( inppar->tasknum == SURFDENSPROF ) && ( inppar->output ) ) {

        real factor, partdens, norm;

        if ( inppar->periodic ) {
            real slicevol;
            real fctr;

            slicevol = ONE;
            for ( i=0; i<DIM; i++ )
                if ( inppar->direction != i )
                    slicevol *= inppar->pbc[i];

            partdens = (real) nref / ( inppar->pbc[0] * inppar->pbc[1] * inppar->pbc[2]);
            factor = (real) counter * partdens * slicevol * drdprof;

            // the one below gives almost the same value
            // fctr = (real) counter * nref / (real) ndprof; // * inppar->profileres;
        }
        else {
            factor = (real) nref * (real) counter;
        }

        printf("particle density:     %21.10f\n", partdens);
        printf("normalization factor: %21.10f\n", factor);

        FILE *fdprof;
        char tmp[MAXSTRLEN];

        sprintf(tmp, "%s%s", inppar->outputprefix, "densprof.dat");
        fdprof = fopen(&tmp[0], "w");

        int hndprof = ndprof / 2;
        for ( i=-hndprof; i<hndprof; i++ ) {

            norm = factor;

            fprintf ( fdprof, "%21.10f %21.10f %21.10f\n", BOHR*i*drdprof, densprof[i+hndprof], densprof[i+hndprof] / norm);
        }

        fclose ( fdprof );

    }

    if ( inppar->tasknum == SURFDENSPROF )
        free ( densprof );

    free(mask);
    free(refmask);
    free(atoms);

    if ( ! ( inppar->xdrread ) )
        fclose ( fxmol );
#ifdef XDRCAP
    else if ( inppar->xdrread )
        xdrfile_close ( xd_read );
#endif

    return 0;
}
