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
#ifdef HAVE_XDRFILE
#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>
#endif
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
#ifdef HAVE_XDRFILE
    XDRFILE * xd_read;
#endif
    FILE * fxmol;
    // FILE * frepxyz;
    FILE * adata;
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
    int adatlnlen;
    int adddata = 0;
    real ntotarea = ZERO;
    real ntotvol = ZERO;
    real * rathist;

    int * mask;
    int * refmask;

    real *zetatom;

    atom_t * atoms;

    i = 0;
    nref = 0;
    nmask = 0;

    if ( strstr(inppar->addeddata, EMPTY) == NULL )
        adddata = 1;
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

#ifdef HAVE_XDRFILE
        xd_read = xdrfile_open ( &inppar->trajectory[0], "r" );
        result_xtc = read_xtc_natoms( &inppar->trajectory[0], &natoms_xtc);

        if ( exdrOK != result_xtc ) {
            printf("Something went wrong opening XDR file\n"); // Error
            exit ( 1 );
        }

        if ( natoms_xtc != natoms ) {
            printf("XDR file and xyz file do not contain the same number of atoms\nNot continuing\n");
            exit ( 1 );
        }
#endif

        fclose ( fxmol );

    }
    else {
        fxmol = fopen(&inppar->trajectory[0], "r");
        natoms = atoi( fgets ( text, MAXSTRLEN, fxmol ) );
        snapsize = xmol_snap_bytesize(fxmol);
        read_xmol(inppar->trajectory, &atoms);
    }

    if ( adddata ) {
        adata = fopen(&inppar->addeddata[0], "r");
        adatlnlen = strlen ( fgets ( text, MAXSTRLEN, adata ) );
    }

    nmask = get_mask(&(mask), inppar->maskkind, inppar->mask, inppar->nkinds, atoms, natoms);

    if  (inppar->tasknum == SURFDENSPROF ) {
            if ( ( strstr ( inppar->refmask, EMPTY ) != NULL ) && ( inppar->nofrags ) )  {
                print_error ( MISSING_INPUT_PARAM, "refmask or fragments" );
                exit ( MISSING_INPUT_PARAM );
            } else if ( ( strstr ( inppar->refmask, EMPTY ) == NULL ) && ( ! ( inppar->nofrags ) ) )  {
                print_error ( CONFLICTING_OPTIONS, "refmask or fragments" );
                exit ( CONFLICTING_OPTIONS );
            }
    }

    int * frags[inppar->numfrags];
    int * frag;
    int ntotfrag = 0;
    int o;
    char buff[10] = "indices";
    int tnt = 0;

    if ( inppar->nofrags ) {
        printf("Using indices given in 'refmask'\n");
        nref = get_mask(&refmask, inppar->refmaskkind, inppar->refmask, inppar->refnkinds, atoms, natoms);
        tnt = nref;

        if ( !(nref) ) {
            printf("Cannot continue with 0 reference atoms\n");
            exit ( 1 );
        }
    }
    else {

        /* check here and wrap fragments to central box */

        printf("Using fragments given in 'fragments'\n");

        for ( o=0; o<inppar->numfrags; o++ ) {
            // int tmpfrg  = get_mask(&frag, &(buff[0]), inppar->fragments[o], inppar->natomsfrag[o], atoms, natoms);
            // check here, and change this later, because it might get ugly with memory freeing
            frags[o] = inppar->fragments[o];
            ntotfrag += inppar->natomsfrag[o];
        }
        tnt = inppar->numfrags;
    }

    int a;
    zetatom = malloc ( natoms * sizeof ( real ) );

    if ( inppar->zetalloc ) {
        for ( a=0; a<natoms; a++ )
            zetatom[a] = inppar->zeta[atoms[a].number];
    }
    else {
        for ( a=0; a<natoms; a++ )
            zetatom[a] = inppar->zetadef;
    }

    if ( inppar->masslloc )
        for ( a=0; a<natoms; a++ )
            atoms[a].mass = inppar->mass[atoms[a].number];

    real *densprof;
    real mxdim;
    int ndprof;
    real drdprof = inppar->profileres;
    real dv;

    if ( inppar->tasknum == SURFDENSPROF ) {

        int mini;

        if ( ! ( inppar->pbcset ) ) {
             print_error ( MISSING_INPUT_PARAM, "pbc" );
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

        if ( adddata )
            rathist = (real *) calloc (ndprof, sizeof(real));
    }

    int seekpoint;
    int ncol = inppar->adatacolstop - inppar->adatacolstrt;
    real adatarray[tnt][ncol];

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
    real *dx;
#ifdef HAVE_XDRFILE
    matrix box;
#endif
    real pbc[DIM];

    for ( i=0; i<DIM; i++ )
        pbc[i] = inppar->pbc[i];

    for ( i=inppar->start; i<inppar->stop; i += inppar->stride )
    {
        /* read snapshot */
        if ( inppar->xdrread ) {
#ifdef HAVE_XDRFILE
            read_xtr_forward ( xd_read, frwrd, atoms, natoms, &box );
            frwrd = inppar->stride;

            int k;
            for ( k=0; k<DIM; k++ )
                pbc[k] = box[k][k];
#endif
        }
        else {
            xmolreader(fxmol, snapsize, i, atoms, natoms);
        }

        if ( adddata ) 
        {

            int nfrg;
            if ( inppar->nofrags )
                nfrg = nref;
            else
                nfrg = inppar->numfrags;

            // if ( strstr ( inppar->task, "rdf" ) != NULL ) {
            //     printf("Sorry, RDFs plus added data currently not implemented\n");
            //     exit ( 1 );
            // }

            seekpoint = counter * adatlnlen * nfrg;
            read_adata(adata, inppar->adatacolstrt, inppar->adatacolstop, &adatarray[0][0], nfrg, seekpoint);
        }

        if ( ( inppar->tasknum == SURFDIST ) || ( inppar->tasknum == SURFDENSPROF ) ) {

            // FU| check here, still need surface area determination

            cube_t surface;
            int fake_n[DIM];
            real fake_origin[DIM];
            real fake_boxv[DIM][DIM];
            real disthi, distlo;
            int * surf_inds;
            real vol;
            int inthi, intlo;
            char tmp[MAXSTRLEN];

            char opref[MAXSTRLEN];
            sprintf(opref, "%s%i_", inppar->outputprefix, i);

            int * direction;
            real * grad;
            real * dstnc;
            real ** surfpts;
            int newsurf= 0;
            int nsurf = 0;

            if ( inppar->load_surface ) {

                FILE *fsurf;

                sprintf(tmp, "%s%i_%s", inppar->loadprefix, i, "atrep_surface.xyz");
                fsurf = fopen(&tmp[0], "r");

                real area;

                if((fgets(text, MAXSTRLEN, fsurf)) != NULL)
                {
                    nsurf = atoi(text);
                    surfpts = (real **) malloc ( nsurf * sizeof ( real *) );

                    int g;
                    for ( g=0; g<nsurf; g++ )
                        surfpts[g] = ( real * ) malloc ( DIM * sizeof ( real ) );

                    direction =  ( int * ) malloc ( nsurf * sizeof ( int ) );
                    grad = ( real * ) malloc ( nsurf * sizeof ( real ) );
                }

                if((fgets(text, MAXSTRLEN, fsurf)) != NULL) {
                    char *dummy;

                    dummy = strtok ( text, " \n" );
                    vol = atof ( dummy );
                    ntotvol += vol;

                    dummy = strtok ( NULL, " \n" );
                    area = atof ( dummy );
                    ntotarea += area;
                }

                int g;
                char *dummy;
                for ( g=0; g<nsurf; g++ )
                    if((fgets(text, MAXSTRLEN, fsurf)) != NULL) {

                        dummy = strtok ( text, " \n" );

                        int k;
                        for ( k=0; k<DIM; k++ ) {
                            dummy = strtok ( NULL, " \n");
                            surfpts[g][k] = atof ( dummy ) / BOHR;
                        }

                        dummy = strtok ( NULL, " \n");
                        grad[g] = atof ( dummy );

                        dummy = strtok ( NULL, " \n");
                        direction[g] = atoi ( dummy );
                    }

                fclose(fsurf);

#ifdef DEBUG
                if ( inppar->surfxyz ) {
                    FILE *fsxyzal;

                    sprintf(tmp, "%s%i_%s", inppar->outputprefix, i, "atrep_surface_rewrite.xyz");
                    fsxyzal = fopen(&tmp[0], "w");;

                    int a, g, k;

                    fprintf ( fsxyzal, "%i\n", nsurf );
                    fprintf ( fsxyzal, "%14.8f %14.8f\n", vol, area );

                    for ( g=0; g<nsurf; g++ ) {
                        fprintf(fsxyzal, "%5s", "X");
                        for ( k=0; k<DIM; k++ ) {
                            fprintf ( fsxyzal, "    %21.10f", BOHR * surfpts[g][k]);
                        }
                        fprintf( fsxyzal, "    %21.10f    %5i\n", grad[g], direction[g]);
                        //fprintf( fsxyzal, "\n" );
                    }

                    fclose ( fsxyzal );
                }
#endif

            }
            else {
                surface = instant_surface_periodic ( mask, atoms, natoms, zetatom, inppar->surfacecutoff, inppar->output, opref, pbc, inppar->resolution, inppar->accuracy, 0, fake_origin, fake_n, fake_boxv, inppar->periodic, 0 );

                if ( inppar->postinterpolate > 1 ) {
                    if ( ! ( inppar->localsurfint ) ) {
                        cube_t fine;

                        if ( inppar->interpolkind == INTERPOLATE_TRILINEAR )
                            fine = interpolate_cube_trilinear ( &surface, inppar->postinterpolate, inppar->periodic );
#ifdef HAVE_EINSPLINE
                        else if ( inppar->interpolkind == INTERPOLATE_BSPLINES )
                            fine = interpolate_cube_bsplines ( &surface, inppar->postinterpolate, inppar->periodic );
#endif

                        free ( surface.atoms );
                        free ( surface.voxels );

                        surface = fine;

                        surface.atoms = fine.atoms;
                        surface.voxels = fine.voxels;

                        if ( inppar->output > 2 ) {
                            sprintf(tmp, "%s%i_%s", inppar->outputprefix, i, "interpolated-instant-surface.cube");
                            write_cubefile(tmp, &surface);
                        }
                    }
                }

                if ( ( inppar->normalization == NORM_BULK) || ( inppar->normalization == NORM_SLAB ) ) {
                    vol = get_bulk_volume ( &surface, inppar->surfacecutoff );

                    if ( inppar->normalization == NORM_SLAB )
                        vol = pbc[0]*pbc[1]*pbc[2] - vol;
                }
                else
                    vol = pbc[0]*pbc[1]*pbc[2];

                // check here, depending on whether we want to look at stuff in the non-solvent phase or in the solvent phase we need to take different volumes
                // printf("%21.10f%21.10f\n", vol, pbc[0]*pbc[1]*pbc[2]);

                ntotvol += vol;

                /* check here, and remove hard-coded surface direction */

                real area = ZERO;

                surfpts = get_2d_representation_ils ( &nsurf, &direction, &grad, &surface, inppar->surfacecutoff, newsurf, surf_inds, inppar->direction, &area, inppar->periodic );

                ntotarea += area;

                // use function write_combined_xmol
                if ( inppar->surfxyz ) {
                    FILE *fsxyzal;

                    sprintf(tmp, "%s%i_%s", inppar->outputprefix, i, "atrep_surface.xyz");
                    fsxyzal = fopen(&tmp[0], "w");;

                    int a, g, k;

                    fprintf ( fsxyzal, "%i\n", nsurf );
                    fprintf ( fsxyzal, "%14.8f %14.8f\n", vol, area );

                    for ( g=0; g<nsurf; g++ ) {
                        fprintf(fsxyzal, "%5s", "X");
                        for ( k=0; k<DIM; k++ ) {
                            fprintf ( fsxyzal, "    %21.10f", BOHR * surfpts[g][k]);
                        }
                        fprintf( fsxyzal, "    %21.10f    %5i\n", grad[g], direction[g]);
                        //fprintf( fsxyzal, "\n" );
                    }

                    fclose ( fsxyzal );
                }
            }

            if ( inppar->tasknum == SURFDENSPROF ) {

                int r;
                signed int ind;
                int hndprof = ndprof / 2.;
                int *fakemask;
                int fakenum;

                if ( counter == 0 )
                    dx = get_box_volels ( &surface );

                int nfrg;

                if ( inppar->nofrags )
                    nfrg = nref;
                else
                    nfrg = inppar->numfrags;

                dstnc = ( real * ) malloc ( nfrg * sizeof ( real ) );

#ifdef OPENMP
                // each thread should have about the same amount of work, so the atomic update will not be too harmful
#pragma omp parallel for default(none) \
                private(r,fakemask,fakenum,ind) \
                shared(dstnc,nfrg,refmask,frags,inppar,surface,direction,natoms,opref,densprof,hndprof,atoms,nsurf,grad,surfpts,pbc,ncol,rathist,adatarray,adddata)
#endif
                for ( r=0; r<nfrg; r++ ) {

                    if ( inppar->nofrags ) {
                        fakemask = &(refmask[r]);
                        fakenum = 1;
                    }
                    else {
                        fakemask = frags[r];
                        fakenum = inppar->natomsfrag[r];
                    }

                    int mnnd;
                    dstnc[r] = get_distance_to_surface ( &mnnd, nsurf, surfpts, direction, grad, atoms, fakemask, fakenum, natoms, pbc, inppar->output, opref, inppar->surfacecutoff, inppar->periodic );
                    // dstnc[r] = get_distance_to_surface ( &mnnd, &surface, nsurf, surfpts, direction, grad, atoms, fakemask, fakenum, natoms, pbc, inppar->output, opref, inppar->surfacecutoff, inppar->periodic );

                    real clspt[DIM];
                    int g;

                    for ( g=0; g<DIM; g++)
                        clspt[g] = surfpts[mnnd][g];

                    if ( ( inppar->postinterpolate > 1 ) && ( inppar->localsurfint ) && ( fabs ( dstnc[r] ) < inppar->ldst ) ) {


                        // need index of point on surface
                        // this is mnnd

                        cube_t fine = local_interpolation ( &surface, surfpts[mnnd], inppar->lint, inppar->interpolkind, inppar->postinterpolate, inppar->outputprefix, pbc, inppar->periodic );
                        // call local interpolation routine

                        int * drctn;
                        real * grd;
                        int nwsrf = 0;
                        int nsrf = 0;
                        int * srf_nds;

                        real rea = ZERO;
                        real ** srfpts;

                        srfpts = get_2d_representation_ils ( &nsrf, &drctn, &grd, &fine, inppar->surfacecutoff, nwsrf, srf_nds, inppar->direction, &rea, 0 );

#ifdef DEBUG
                        if ( inppar->surfxyz ) {
                            FILE *fsxyzal;

                            sprintf(tmp, "%s%i_%s", inppar->outputprefix, r, "test_atrep_surface.xyz");
                            fsxyzal = fopen(&tmp[0], "w");;

                            int a, g, k;

                            fprintf ( fsxyzal, "%i\n\n", nsrf );

                            for ( g=0; g<nsrf; g++ ) {
                                fprintf(fsxyzal, "%5s", "X");
                                for ( k=0; k<DIM; k++ ) {
                                    fprintf ( fsxyzal, "    %21.10f", BOHR * srfpts[g][k]);
                                }
                                fprintf( fsxyzal, "\n" );
                            }

                            fclose ( fsxyzal );
                        }
#endif
                        // get distance again

                        int mnnd;

                        dstnc[r] = get_distance_to_surface ( &mnnd, nsrf, srfpts, drctn, grd, atoms, fakemask, fakenum, natoms, pbc, inppar->output, opref, inppar->surfacecutoff, inppar->periodic );
                        // dstnc[r] = get_distance_to_surface ( &mnnd, &fine, nsrf, srfpts, drctn, grd, atoms, fakemask, fakenum, natoms, pbc, inppar->output, opref, inppar->surfacecutoff, inppar->periodic );

                        for ( g=0; g<DIM; g++)
                            clspt[g] = surfpts[mnnd][g];

                        int l;
                        for ( l=0; l<nsrf; l++ )
                            free ( srfpts[l] );

                        if ( nwsrf )
                            free ( srf_nds );

                        free ( grd );
                        free ( srfpts );
                        free ( drctn );

                        free ( fine.atoms );
                        free ( fine.voxels );
                    }

                    ind = ( int ) floor ( dstnc[r] / inppar->profileres );

                    real tmpdat;

                    //FUDO| should we do other stuff here as well, more options for ncol?
                    if ( adddata ) {
                        if ( ncol == 3 ) {
                            real dstncvec[DIM];
                            real tmpcom[DIM];

                            get_center_of_mass ( tmpcom, atoms, fakemask, fakenum);

                            real nrm_a = 0.;
                            real nrm_b = 0.;

                            if ( inppar->periodic )
                                get_distance_vector_periodic ( dstncvec, clspt, tmpcom, pbc );
                            else
                                for ( g=0; g<ncol; g++ )
                                    dstncvec[g] = clspt[g] - tmpcom[g];

                            // FUDO| periodicity 
                            for ( g=0; g<ncol; g++ ) {
                               // dstncvec[g] = clspt[g] - tmpcom[g];
                               nrm_a += dstncvec[g]*dstncvec[g];
                               nrm_b += adatarray[r][g]*adatarray[r][g];
                            }

                            nrm_a = sqrt(nrm_a);
                            if ( fabsf ( nrm_a - dstnc[r] ) > 1.e-06 )
                                printf("Distances not congruent, something's wrong %14.8f %14.8f\n", nrm_a, dstnc[r]);

                            nrm_b = sqrt(nrm_b);

                            // printf("%5i %14.8f %14.8f\n", r, dstncvec[0], adatarray[r][0]);
                            // printf("%5i %14.8f %14.8f\n", r, dstncvec[1], adatarray[r][1]);
                            // printf("%5i %14.8f %14.8f\n", r, dstncvec[2], adatarray[r][2]);

                            tmpdat = dstncvec[0] * adatarray[r][0] + dstncvec[1] * adatarray[r][0] + dstncvec[2] * adatarray[r][2];
                            tmpdat /= (nrm_a * nrm_b);
                        }
                        else {
                            tmpdat = adatarray[r][0];
                            // printf("%i %14.8f\n", r, adatarray[r][0]);
                        }
                    }

#pragma omp atomic update
                    densprof[ hndprof + ind ] += 1.; // / nsurf;
                    if ( adddata )
#pragma omp atomic update
                        rathist [ hndprof + ind ] += tmpdat;

                }

                if ( inppar->output ) {
                    sprintf(tmp, "%s%i_%s", inppar->outputprefix, i, "surfdist.dat");
                    fsdist = fopen(&tmp[0], "w");

                    fprintf ( fsdist, "#            index              distance [a.u.]\n");
                    for ( r=0; r<nfrg; r++ )
                        fprintf ( fsdist, "%21i %21.10f\n", r, dstnc[r]);

                    fclose ( fsdist );
                }

            }

            /* check here, and move stuff for refinement box creation somewhere else */

            int k;
            for ( k=0; k<nsurf; k++ )
                free ( surfpts[k] );

            if ( newsurf )
                free ( surf_inds );

            free ( grad );
            free ( dstnc );
            free ( surfpts );
            free ( direction );
            if ( !( inppar->load_surface ) ) {
                free ( surface.atoms );
                free ( surface.voxels );
            }
        }

        printf("%4.2f %% done\r", (real) counter / ntot * 100.);
        fflush(stdout);

        counter++;

    }
    printf("%4.2f %% done\n", (real) counter / ntot * 100.);

    if ( ( inppar->tasknum == SURFDENSPROF ) && ( inppar->output ) ) {

        int natdens;
        int navsurf;
        real avarea;
        real avvol;
        real factor, partdens, norm;

        if ( inppar->nofrags )
            natdens = nref;
        else
            natdens = inppar->numfrags;

        real smarea;
        if ( inppar->periodic ) {
            real fctr;

            avarea = ntotarea / counter;
            avvol = ntotvol / counter;
            smarea = avarea;

            partdens = (real) natdens / ( avvol );
            // partdens = (real) natdens / ( pbc[0] * pbc[1] * pbc[2]);

            factor = (real) counter * drdprof * smarea * partdens;

        }
        else {
            factor = (real) counter * drdprof * smarea * natdens;
        }

        printf("particle density:     %21.10f\n", partdens);
        printf("normalization factor: %21.10f\n", factor);
        printf("average surface area: %21.10f\n", smarea);

        FILE *fdprof;
        FILE *fadatout;

        char tmp[MAXSTRLEN];

        sprintf(tmp, "%s%s", inppar->outputprefix, "densprof.dat");
        fdprof = fopen(&tmp[0], "w");

        if ( adddata ) {
            sprintf(tmp, "%s%s", inppar->outputprefix, "added.dat");
            fadatout = fopen(&tmp[0], "w");
        }

        int hndprof = ndprof / 2;
        real hdrdprof = drdprof / 2.;
        for ( i=-hndprof; i<hndprof; i++ ) {

            norm = factor;

            fprintf ( fdprof, "%21.10f %21.10f %21.10f\n", BOHR*i*drdprof+hdrdprof, densprof[i+hndprof], densprof[i+hndprof] / norm);
            if ( adddata )
                if ( densprof[i+hndprof] / norm > 1.e-8 )
                    fprintf ( fadatout, "%21.10f %21.10f %21.10f\n", BOHR*i*drdprof+hdrdprof, rathist[i+hndprof], rathist[i+hndprof] / densprof[i+hndprof]);
                else
                    fprintf ( fadatout, "%21.10f %21.10f %21.10f\n", BOHR*i*drdprof+hdrdprof, 0., 0.);
        }

        fclose ( fdprof );
        if ( adddata )
            fclose ( fadatout );

    }

    if ( ! ( inppar->nofrags ) ) {
        for ( i=0; i<inppar->numfrags; i++ ) {
            // free(inppar->fragments[i]);
            free(frags[i]);
        }

        free(inppar->natomsfrag);
        free(inppar->fragments);
    }
    else {
        free ( refmask );
    }

    if ( inppar->tasknum == SURFDENSPROF ) {
        free ( densprof );
        free ( dx );
    }

    free ( zetatom );
    if ( inppar->zetalloc ) {
        free ( inppar->zeta );
    }

    if ( inppar->masslloc )
        free ( inppar->mass );

    if ( adddata )
        free(rathist);

    free(mask);
    // free(refmask);
    free(atoms);

    if ( ! ( inppar->xdrread ) )
        fclose ( fxmol );
#ifdef HAVE_XDRFILE
    else if ( inppar->xdrread )
        xdrfile_close ( xd_read );
#endif

    return 0;
}

// void read_adata(FILE *adata, int colstrt, int colstop, real **adatarray, int nlines, int seekpoint)
void read_adata(FILE *adata, int colstrt, int colstop, real *adatarray, int nlines, int seekpoint)
{
    int i, j;
    char text[MAXSTRLEN];
    char * dummy;
    int col;

    int ncol = colstop - colstrt;
    int baseind = 0;

    fseek ( adata, seekpoint, SEEK_SET );

    for ( i=0; i<nlines; i++ )
    {
        /* get line */ 

        fgets ( text, MAXSTRLEN, adata );

        // printf("%s\n", text);
        // adatarray[i] = atof(text);

        /* split line adatacol times */ 

        dummy = strtok(text, " \n");

        for ( j=0; j<colstrt; j++ )
            dummy = strtok ( NULL, " \n");

        //FUDO| maybe just use the two-dimensional array?

        int cnt = 0;
        for ( j=colstrt; j<colstop; j++ )
        {
            adatarray[baseind+cnt] = atof(dummy);
            // adatarray[i][cnt] = atof(dummy);
            dummy = strtok ( NULL, " \n");
            cnt++;
        }

        // printf("\n");
        // printf("read: %s\n", dummy);

        /* copy data to array */

        // adatarray[i] = atof(dummy);
        
        baseind += ncol;
    }

}
