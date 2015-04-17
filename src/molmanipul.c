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

void distance_matrix(atom_t * atoms, int * mask, int natoms, real ** distmat, int periodic, real *pbc)
{
    int i, j;
    real dist;

    for ( i=0; i<natoms; i++ )
    {
        for ( j=i; j<natoms; j++ )
        {
        // calculate distance also between the same atoms, no significant overhead, IMHumbleO
        // but this is not periodic, damn, but it'd better be...
        if ( periodic )
            dist = get_distance_periodic ( &(atoms[mask[i]].coords[0]), &(atoms[mask[j]].coords[0]), pbc );
        else
            dist = get_distance ( &(atoms[mask[i]].coords[0]), &(atoms[mask[j]].coords[0]) );

        distmat[i][j] = dist;
        distmat[j][i] = dist;

        }
    }
}

// check here and include mask in bond search
int ** guess_fragments(atom_t * atoms, int natoms, int * mask, int * numfrags, real blfudge, int periodic, real *pbc)
{
    int i, j;
    int * countbonds;
    char ** frags;
    int ** bondinds;

    real accbond;
    real ** distmat;
    real ** acceptable_bondlengths;

    char ** bonds;
    char dummy[MAXSTRLEN];

    // real acceptable_bondlengths[LAST_ATOM][LAST_ATOM];

    acceptable_bondlengths = ( real ** ) malloc ( LAST_ATOM * sizeof ( real * ) );

    for ( i=0; i<LAST_ATOM; i++ )
    {
        acceptable_bondlengths[i] = ( real * ) malloc ( LAST_ATOM * sizeof ( acceptable_bondlengths[i]) );
    }

    countbonds = ( int * ) malloc ( natoms * sizeof ( int ) );
    distmat = (real **) malloc ( natoms * sizeof ( real *) );
    bonds = ( char ** ) malloc ( natoms * sizeof ( char * ) );
    bondinds = ( int ** ) malloc ( natoms * sizeof ( int * ) );
    frags = ( char ** ) malloc ( natoms * sizeof ( char * ) );

    int * isprocessed = (int *) malloc ( natoms * sizeof(int) );

    for ( i=0; i<natoms; i++ )
    {
        distmat[i] = (real *) malloc ( natoms * sizeof ( distmat[i] ) );
        bonds[i] = ( char * ) malloc ( MAXSTRLEN * sizeof ( bonds[i] ) );
        frags[i] = ( char * ) malloc ( MAXSTRLEN * sizeof ( bonds[i] ) );

        countbonds[i] = 0;
        isprocessed[i] = 0;
        // distmat[i] = (real *) malloc ( natoms * sizeof ( real ) );
    }

    distance_matrix(atoms, mask, natoms, distmat, periodic, pbc);
    guess_acceptable_bond_lengths( acceptable_bondlengths );

    for ( i=0; i<natoms; i++ )
    {
        // check here, we only want the bonds and identify the connected atom by array index
        // sprintf(bonds[i], "%i,", i);
        for ( j=0; j<natoms; j++ )
        // for ( j=i+1; j<natoms; j++ )
        {
            if ( i == j) 
                continue;

            accbond = blfudge * acceptable_bondlengths[atoms[mask[i]].number][atoms[mask[j]].number];

            if ( distmat[i][j] < accbond )
            {
                sprintf(dummy, "%i,", j);
                strcat(bonds[i], dummy);
                countbonds[i]++;
                // printf("found a bond between %i %i\n", i, j);
            }
        }

        // allocate integer arrays with bonded atoms indicees and then concatenate all

        char symdex[8] = "indices";
        bondinds[i] = get_mask(&symdex[0], bonds[i], countbonds[i], atoms, natoms);

        // printf("%s\n", bonds[i]);
        // for ( j=0; j<countbonds[i]; j++ )
        //     printf("%i  ", bondinds[i][j]);
        // printf("\n");
    }

    int ** fragatoms;
    fragatoms = ( int ** ) malloc ( natoms * sizeof ( int * ) );

    *numfrags = 0;

    int fragcount = 0;
    for ( i=0; i<natoms; i++ )
    {
        if ( isprocessed[i] )
            continue;

        (*numfrags)++;

        int natomsfrag = 1;

        char fragmask[MAXSTRLEN];

        sprintf( fragmask, "%i,", i );

        recurse_concatenate ( &fragmask[0], bondinds, countbonds, i, isprocessed, &natomsfrag);

        char symdex[8] = "indices";
        fragatoms[fragcount] = get_mask( &symdex[0], fragmask, natomsfrag, atoms, natoms);

        bubble_sort_ints(fragatoms[fragcount], natomsfrag);

        fragcount++;
    }

#ifdef DEBUG
    for ( i=0; i<natoms; i++ )
    {
        for ( j=0; j<natoms; j++ )
        {
            printf("%f  ", distmat[i][j]);
        }
        printf("\n");
    }
#endif

    // check here, and deallocate all the arrays
    // and also return array of arrays of indices instead of character mask

    for ( i=0; i<LAST_ATOM; i++ )
    {
        free ( acceptable_bondlengths[i] );
    }

    free(acceptable_bondlengths);

    for ( i=0; i<natoms; i++ )
    {
        free(distmat[i]);
        free(bonds[i]);
        free(frags[i]);
        free(bondinds[i]);
    }

    free(countbonds);
    free(distmat);
    free(bonds);
    free(bondinds);
    free(frags);
    free(isprocessed);

    return fragatoms;
}

void recurse_concatenate ( char * frags, int ** bondinds, int * countbonds, int atomind, int * isprocessed, int * natomsfrag)
{
    int i;

    for ( i=0; i<countbonds[atomind]; i++ )
    {
        if ( isprocessed[bondinds[atomind][i]] )
        // if ( ( isprocessed[bondinds[atomind][i]] ) || ( isprocessed[atomind] ))
        {
            continue;
        }

        char dummy[MAXSTRLEN];

        // can this by accident flooded with previous data?
        sprintf(dummy, "%i,", bondinds[atomind][i]);
        strcat ( frags, dummy );

        (*natomsfrag)++;
        isprocessed[atomind] = 1;

        recurse_concatenate ( frags, bondinds, countbonds, bondinds[atomind][i], isprocessed, natomsfrag );

        isprocessed[bondinds[atomind][i]] = 1;
    }
}

void get_center_of_mass ( real * com, atom_t * atoms, int * mask, int natoms )
{
    // real *com;
    real totmass;
    int i, j;

    // com = (real *) calloc ( DIM, sizeof ( real ) );
    for ( i=0; i<DIM; i++ )
        com[i] = ZERO;

    totmass = ZERO;

    for ( i=0; i<natoms; i++ ) {
        for ( j=0; j<DIM; j++ )
            com[j] += atoms[mask[i]].mass * atoms[mask[i]].coords[j];

        totmass += atoms[mask[i]].mass;
    }

    for ( i=0; i<DIM; i++ )
        com[i] /= totmass;

    // return com;
}

void move_atoms ( atom_t * atoms, int * mask, int natoms, real * mover )
{
    int i, j;

    for ( i=0; i<natoms; i++ )
        for ( j=0; j<DIM; j++ )
            atoms[mask[i]].coords[j] += mover[j];

}
