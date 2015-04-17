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
#include "cube.h"
#include "atom_param.h"
#include "io.h"
#include "errors.h"
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <math.h>
#include <zlib.h>

void parse_input_file(input_t * inppar, char * inpfile)
{
    FILE *datei;

    if((datei = fopen(inpfile, "r")))
    {
        char text[MAXSTRLEN];
        char * variable;
        char * value;

        printf("\nParsing input file.\n");
        printf("\nRepetition of input parameters:\n-----------------------------------------\n");
        while(fgets(text, MAXSTRLEN, datei) != NULL)
        {
            if (strstr(text, "=") != NULL)
            {
            variable = strtok (text, "=\n");
            value = strtok (NULL, "=\n");

            // something like that, but needs changes in strtok (because now the first character would not necessarily be #)
            // if (strncmp(variable, "#", 1) == 0)
            // {
            //     break;
            // }

            set_input_value(inppar, variable, value);
            }
        }
        fclose(datei);
        printf("-----------------------------------------\n\n");
    }
    else
    {
        print_error(FILE_NOT_FOUND, inpfile);
        exit(1);
        return;
    }
}

void set_input_value(input_t *inppar, char *variable, char *value)
{

    long key;
    char * dummy;

    // check here and exclude lines that start with a hash
    if ( strncmp ( variable, "#", 1 ) == 0 )
        return;
            
    dummy = strtok(variable, " ");
    variable = dummy;
    int counter = 0;

    for (key = 0; key < NUMKEYS; key++)
    {
        if (strncmp(variable, keywords[key], 20) == 0)
        {
            counter++;
            dummy = strtok(value, " ");
            value = dummy;

            printf("%s is set to: '%s'\n", variable, value);
            if ((keywords[key] == "task"))
            {
                strcpy(inppar->task, value);
                inppar-> tasknum = assign_task ( inppar->task );
            }
            else if ((keywords[key] == "structure"))
            {
                strcpy(inppar->structure, value);
            }
            else if ((keywords[key] == "oprefix"))
            {
                strcpy(inppar->outputprefix, value);
            }
            else if ((keywords[key] == "surfrefinement"))
            {
                inppar->surfrefinement = atof(value);
            }
            else if ((keywords[key] == "refinementitpl"))
            {
                inppar->refineitpl = atoi(value);
            }
            else if ((keywords[key] == "surfacecutoff"))
            {
                inppar->surfacecutoff = atof(value);
            }
            else if ((keywords[key] == "zeta"))
            {
                inppar->zeta = atof(value);
            }
            else if ((keywords[key] == "postinterpolate"))
            {
                inppar->postinterpolate = atoi(value);
            }
            else if ((keywords[key] == "roughsurf"))
            {
                inppar->roughsurf = atoi(value);
            }
            else if ((keywords[key] == "pbc"))
            {
                inppar->pbcset = 1;
                /* check here, maybe we don't always want to set inppar->periodic, when we set pbc, but I dunno */
                // inppar->periodic = 1;
                int i=0;
                real conv = 1.;

                if ( strstr ( dummy, "ANG" ) != NULL ) {
                    conv = ANG2BOHR;
                    dummy = strtok ( NULL, " " );
                }
                else if ( strstr ( dummy, "NM" ) != NULL ) {
                    conv = NM2BOHR;
                    dummy = strtok ( NULL, " " );
                }

                while ( dummy != NULL )
                {
                    inppar->pbc[i] = atof(dummy);
                    inppar->pbc[i] *= conv;

                    dummy = strtok(NULL, " ");
                    i++;
                }

                if ( i == 1 )
                    for ( i=1; i<DIM; i++ )
                        inppar->pbc[i] = inppar->pbc[0];

            }
            else if ((keywords[key] == "solcenter"))
            {
                inppar->solset = 1;
                int i=0;

                while ( dummy != NULL )
                {
                    inppar->solcenter[i] = atof(dummy);

                    dummy = strtok(NULL, " ");
                    i++;
                }

                if ( i == 1 )
                    for ( i=1; i<DIM; i++ )
                        inppar->solcenter[i] = inppar->solcenter[0];

            }
            else if ((keywords[key] == "refcenter"))
            {
                int i=0;

                while ( dummy != NULL )
                {
                    inppar->refcenter[i] = atof(dummy);

                    dummy = strtok(NULL, " ");
                    i++;
                }

                inppar->refcenterset = 1;
            }
            else if ((keywords[key] == "mask"))
            {
                strcpy(inppar->maskkind, value);

                dummy = strtok(NULL, " ");
                inppar->nkinds = 0;

                value = dummy;
                strcpy(inppar->mask, value);

                while ( dummy != NULL )
                {
                    dummy = strtok(NULL, " ");
                    if ( dummy != NULL)
                    {
                        strcat(inppar->mask, ",");
                        strcat(inppar->mask, dummy);
                    }
                    inppar->nkinds++;
                }
            }
            else if ((keywords[key] == "refmask"))
            {
                strcpy(inppar->refmaskkind, value);

                dummy = strtok(NULL, " ");
                inppar->refnkinds = 0;

                value = dummy;
                strcpy(inppar->refmask, value);

                while ( dummy != NULL )
                {
                    dummy = strtok(NULL, " ");
                    if ( dummy != NULL)
                    {
                        strcat(inppar->refmask, ",");
                        strcat(inppar->refmask, dummy);
                    }
                    inppar->refnkinds++;
                }
            }
            else if ((keywords[key] == "batch"))
            {
                inppar->start = atoi(value);
                inppar->batchmode = 1;

                dummy = strtok(NULL, " ");
                inppar->stop = atoi(dummy);

                dummy = strtok(NULL, " ");
                inppar->stride = atoi(dummy);
            }
            else if ((keywords[key] == "output"))
            {
                if (strstr(value, "silent") != NULL)
                {
                    inppar->output = 0;
                    printf("silent mode activated\n");
                }
                else
                    if (strstr(value, "normal") != NULL)
                    {
                        inppar->output = 1;
                    }
                    else
                        if (strstr(value, "debug") != NULL)
                        {
                            inppar->output = 2;
                        }
            }
            else if ((keywords[key] == "fragments"))
            {
                inppar->nofrags = 0;
                inppar->numfrags = atoi(value);

                inppar->natomsfrag = (int *) malloc(inppar->numfrags * sizeof(int));

                int i;

                for ( i=0; i<inppar->numfrags; i++ )
                {
                    inppar->fragments[i] = (char *) malloc(MAXSTRLEN * sizeof(char));
                }

                i = 0;
                dummy = strtok(NULL, " ");
                inppar->fragtotnatoms = 0;
                while ( dummy != NULL )
                {
                    inppar->natomsfrag[i] = 0;
                    while ( strstr(dummy, ";") == NULL )
                    {
                        strcat((inppar->fragments)[i], dummy);
                        strcat((inppar->fragments)[i], ",");
                        dummy = strtok(NULL, " ");
                        inppar->natomsfrag[i] += 1;
                    }

                    int tmplen = strlen(dummy);
                    char buffer[MAXSTRLEN] = "";

                    strncpy(&buffer[0], dummy, tmplen-1);
                    strcat((inppar->fragments)[i], buffer);
                    inppar->natomsfrag[i] += 1;

                    // put here later, does the last "," hurt the get_mask function?
                    inppar->fragtotnatoms += inppar->natomsfrag[i];
                    
                    printf("%i atoms in fragment #%i: %s\n", inppar->natomsfrag[i], i, inppar->fragments[i]);

                    dummy = strtok(NULL, " ");

                    strcat(inppar->fragmask, (inppar->fragments)[i]);
                    
                    if ( dummy != NULL )
                        strcat(inppar->fragmask, ",");

                    // printf("%s\n", (inppar->fragments)[i]);
                    i++;
                }
            }
            else if ((keywords[key] == "periodic"))
            {
                inppar->periodic = atoi(value);
            }
            else if ((keywords[key] == "wrap"))
            {
                inppar->wrap = atoi(value);
            }
            else if ((keywords[key] == "blfudge"))
            {
                inppar->blfudge = atof(value);
            }
            else if ((keywords[key] == "guessfragments"))
            {
                inppar->guessfragments = atoi(value);
            }
            else if ((keywords[key] == "trajectory"))
            {
                strcpy(inppar->trajectory, value);
                inppar->trajmode = 1;
            }
            else if ((keywords[key] == "resolution"))
            {
                inppar->resolution = atof(value);
            }
            else if ((keywords[key] == "xdrread"))
            {
#ifdef XDRCAP
                inppar->xdrread = atoi(value);
#else
                printf("SURF not compiled with XDR support\n");
                exit ( 1 );
#endif
            }
            else if ((keywords[key] == "ignorefirst"))
            {
                inppar->gnrfrst = atoi(value);
            }
            else if ((keywords[key] == "accuracy"))
            {
                inppar->accuracy = atof ( value );
            }
            else if ((keywords[key] == "direction"))
            {
                inppar->direction = atoi ( value );
            }
            else if ((keywords[key] == "dummy"))
            {
                real conv = 1.;
                int j;

                if ( ! ( inppar->pbcset ) )
                    print_error ( MISSING_INPUT_PARAM, "pbc before dummy" );

                // change that whole shit to a function that figures out conversion factor
                if ( strstr ( dummy, "ANG" ) != NULL ) {
                    conv = ANG2BOHR;
                    dummy = strtok ( NULL, " " );
                }
                else if ( strstr ( dummy, "NM" ) != NULL ) {
                    conv = NM2BOHR;
                    dummy = strtok ( NULL, " " );
                }

                if ( strstr ( dummy, "boxcenter" ) != NULL ) {
                    for ( j=0; j<DIM; j++ )
                        inppar->dumatom.coords[j] = inppar->pbc[j] / 2.;

                    inppar->dummy = BOXCENTER;
                    dummy = strtok ( NULL, " " );
                }

                else if ( strstr ( dummy, "boxedgehigh" ) != NULL ) {
                    for ( j=0; j<DIM; j++ )
                        inppar->dumatom.coords[j] = inppar->pbc[j];

                    inppar->dummy = BOXEDGEHI;
                    dummy = strtok ( NULL, " " );
                }

                else if ( strstr ( dummy, "boxedgelow" ) != NULL ) {
                    for ( j=0; j<DIM; j++ )
                        inppar->dumatom.coords[j] = ZERO;

                    inppar->dummy = BOXEDGELO;
                    dummy = strtok ( NULL, " " );
                }
                else {

                    inppar->dummy = USERDEFINED;
                    int i=0;

                    while ( dummy != NULL )
                    {
                        inppar->dumatom.coords[i] = atof(dummy);
                        inppar->dumatom.coords[i] *= conv;

                        dummy = strtok(NULL, " ");
                        i++;
                    }

                    if ( i == 1 )
                        for ( i=1; i<DIM; i++ )
                            inppar->dumatom.coords[i] = inppar->dumatom.coords[0];
                }

            }
        }

        if ( counter )
            break;
    }

    if ( ! ( counter  ) )
    {
        printf("\nKeyword '%s' is not known to me. Please check your input file again\n\n", variable);
        exit(1);
    }

    return;
}

void write_cubefile(char * filename, cube_t * cube)
{
    FILE * datei;
    int i, j;

    char spacer[]="  ";

    if((datei = fopen(filename, "w")))
    {
        fprintf(datei, "%s", cube->title);
        fprintf(datei, "%s", cube->comment);

        fprintf(datei, "%s%3i", spacer, cube->natoms);
        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%10.6f", spacer, cube->origin[i]);
        }
        fprintf(datei, "\n");
        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%3i", spacer, cube->n[i]);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->boxv[i][j]);
            }
            fprintf(datei, "\n");
        }

        for (i=0; i<cube->natoms; i++)
        {
            fprintf(datei, "%s%3i%s%10.6f", spacer, cube->atoms[i].number, spacer, cube->atoms[i].charge);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->atoms[i].coords[j]);
            }
            fprintf(datei, "\n");
        }

        j = 0;
        for (i=0; i<cube->nvoxels; i++)
        {
            fprintf(datei, "%s%12.5e", spacer, cube->voxels[i].data);
            if (i % 6 == 5)
            {
                fprintf(datei, "\n");
            }
        }
        fprintf(datei, "\n");
    }
    // else
    // {
    //     HELLO
    // }

    fclose(datei);
}

void write_cubefile_indices(char * filename, cube_t * cube, int * indices)
{
    FILE * datei;
    int i, j;

    char spacer[]="  ";

    if((datei = fopen(filename, "w")))
    {
        fprintf(datei, "%s", cube->title);
        fprintf(datei, "%s", cube->comment);

        fprintf(datei, "%s%3i", spacer, cube->natoms);
        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%10.6f", spacer, cube->origin[i]);
        }
        fprintf(datei, "\n");
        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%3i", spacer, cube->n[i]);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->boxv[i][j]);
            }
            fprintf(datei, "\n");
        }

        for (i=0; i<cube->natoms; i++)
        {
            fprintf(datei, "%s%3i%s%10.6f", spacer, cube->atoms[i].number, spacer, cube->atoms[i].charge);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->atoms[i].coords[j]);
            }
            fprintf(datei, "\n");
        }

        j = 0;
        for (i=0; i<cube->nvoxels; i++)
        {
            if (i == indices[j])
            {
                fprintf(datei, "%s%12.5e", spacer, cube->voxels[i].data);
                j++;
            }
            else
            {
                fprintf(datei, "%s%12.5e", spacer, ZERO);
            }

            if (i % 6 == 5)
            {
                fprintf(datei, "\n");
            }
        }
        fprintf(datei, "\n");
    }

    fclose(datei);
}

void write_cubefile_indices_plain(char * filename, cube_t * cube, int * indices)
{
    FILE * datei;
    int i, j;

    char spacer[]="  ";

    if((datei = fopen(filename, "w")))
    {
        fprintf(datei, "%s", cube->title);
        fprintf(datei, "%s", cube->comment);

        fprintf(datei, "%s%3i", spacer, cube->natoms);
        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%10.6f", spacer, cube->origin[i]);
        }
        fprintf(datei, "\n");
        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%3i", spacer, cube->n[i]);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->boxv[i][j]);
            }
            fprintf(datei, "\n");
        }

        for (i=0; i<cube->natoms; i++)
        {
            fprintf(datei, "%s%3i%s%10.6f", spacer, cube->atoms[i].number, spacer, cube->atoms[i].charge);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->atoms[i].coords[j]);
            }
            fprintf(datei, "\n");
        }

        j = 0;
        for (i=0; i<cube->nvoxels; i++)
        {
            if (i == indices[j])
            {
                fprintf(datei, "%s%12.5e", spacer, ONE);
                j++;
            }
            else
            {
                fprintf(datei, "%s%12.5e", spacer, ZERO);
            }
            if (i % 6 == 5)
            {
                fprintf(datei, "\n");
            }
        }
        fprintf(datei, "\n");
    }

    fclose(datei);
}

void write_cubefile_array(char * filename, cube_t * cube, real data[])
{
    FILE * datei;
    int i, j;

    char spacer[]="  ";

    if((datei = fopen(filename, "w")))
    {
        fprintf(datei, "%s", cube->title);
        fprintf(datei, "%s", cube->comment);

        fprintf(datei, "%s%3i", spacer, cube->natoms);
        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%10.6f", spacer, cube->origin[i]);
        }
        fprintf(datei, "\n");
        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%3i", spacer, cube->n[i]);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->boxv[i][j]);
            }
            fprintf(datei, "\n");
        }

        for (i=0; i<cube->natoms; i++)
        {
            fprintf(datei, "%s%3i%s%10.6f", spacer, cube->atoms[i].number, spacer, cube->atoms[i].charge);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->atoms[i].coords[j]);
            }
            fprintf(datei, "\n");
        }

        j = 0;
        for (i=0; i<cube->nvoxels; i++)
        {
            fprintf(datei, "%s%12.5e", spacer, data[i]);
            if (i % 6 == 5)
            {
                fprintf(datei, "\n");
            }
        }
        fprintf(datei, "\n");
    }

    fclose(datei);
}

void set_input_defaults(input_t * inppar)
{
    strcpy(inppar->task, "bader");
    strcpy(inppar->structure, EMPTY);
    strcpy(inppar->outputprefix, "");
    strcpy(inppar->mask, "0,");
    strcpy(inppar->refmask, EMPTY);
    strcpy(inppar->maskkind, "notatoms");
    strcpy(inppar->refmaskkind, "notatoms");
    strcpy(inppar->inpfile, EMPTY);
    inppar->output = 1;
    inppar->surfacecutoff = 0.016* sqr(BOHR) * BOHR;
    inppar->zeta = 2.4/BOHR;
    inppar->roughsurf = 0;
    inppar->batchmode = 0;
    inppar->nkinds = 1;
    inppar->refnkinds = 1;
    inppar->surfrefinement = -50;
    inppar->refineitpl = 0;
    inppar->numfrags = 0;
    inppar->fragtotnatoms = 0;
    inppar->periodic = 0;
    inppar->wrap = 0;
    inppar->blfudge = ONE;
    inppar->guessfragments = 0;
    inppar->nofrags = 1;
    strcpy(inppar->trajectory, EMPTY);
    inppar->trajmode = 0;
    inppar->resolution = 0.1 / BOHR;
    inppar->xdrread = 0;
    inppar->othercenter = 0;
    inppar->rprof_num = 0;
    inppar->accuracy = 1.e-5;
    inppar->tasknum = NOTASK;
    assign_atom_parameters ( "symbol", "Du", &(inppar->dumatom) );
    inppar->dummy = 0;
    inppar->direction = 2;

    int i;
    for ( i=0; i<DIM; i++ )
    {
        inppar->dumatom.coords[i] = ZERO;
        inppar->pbcset = 0;
        inppar->pbc[i] = ZERO;
        inppar->solset = 0;
        inppar->solcenter[i] = ZERO;
        inppar->refcenterset = 0;
        inppar->refcenter[i] = ZERO;
    }
}

int * get_mask(char * maskkind, char * mask, int nkinds, atom_t * atoms, int natoms)
{
    int i, j, k;
    int * indices;
    int * kinds;
    char * dummy;
    int check = 0;

    kinds = (int *) malloc(nkinds * sizeof(int));

    char dummask[MAXSTRLEN];
    strcpy(&(dummask[0]), mask);

    dummy = strtok(dummask, ",");

    for ( i=0; i<nkinds; i++ )
    {
        kinds[i] = atoi(dummy);
        dummy = strtok(NULL, ",");
    }

    if ( strncmp(maskkind, "indices", 5) == 0 )
        indices = (int *) malloc((nkinds+1) * sizeof(int));
    else
        indices = (int *) malloc((natoms+1) * sizeof(int));

    if ( strncmp(maskkind, "atoms", 5) == 0 )
    {
        k = 0;
        for ( i=0; i<natoms; i++)
        {
            for ( j=0; j<nkinds; j++ )
            {
                if ( atoms[i].number == kinds[j] )
                {
                    indices[k] = i;
                    k++;
                    break;
                }
            }
        }
        indices[k] = -1;
    }
    else if ( strncmp(maskkind, "notatoms", 5) == 0 )
    {
        k = 0;
        for ( i=0; i<natoms; i++)
        {
            check = 1;
            for ( j=0; j<nkinds; j++ )
            {
                if ( atoms[i].number == kinds[j] )
                {
                    check = 0;
                    break;
                }
            }
            if ( check )
            {
                indices[k] = i;
                k++;
            }
        }
        indices[k] = -1;
    }
    else if ( strncmp(maskkind, "indices", 7) == 0 )
    {
        for ( i=0; i<nkinds; i++ )
        {
            indices[i] = kinds[i];
        }
        indices[i] = -1;
    }
    else if ( strncmp(maskkind, "notindices", 7) == 0 )
    {
        k = 0;
        for ( i=0; i<natoms; i++ )
        {
            check = 1;
            for ( j=0; j<nkinds; j++ )
            {
                if ( kinds[j] == i ) {
                    check = 0;
                    break;
                }
            }
            if ( check ) {
                indices[k] = i;
                k++;
            }
        }
        indices[k] = -1;
    }

    free(kinds);
    return indices;
}

int read_xmol(char * coordfile, atom_t ** atoms)
{
    FILE * data;
    char text[MAXSTRLEN];
    char * ccoord;
    char * symbol;
    int count=0, i, natoms;

    if((data = fopen(coordfile, "r")))
    {
        if((fgets(text, MAXSTRLEN, data)) != NULL)
        {
            natoms = atoi(text);

            *atoms = (atom_t *)malloc(natoms * sizeof(atom_t));

            count++;
        }
        else
            return -1; 
    }
    else
    {
        printf("Coordinate '%s' file could not be opened\n", coordfile);
        return 0;
    }

    while (fgets(text, MAXSTRLEN, data) != NULL)
    {
        if ((count-2) == natoms)
        {
            break;
        }

        if((count < 2))
        {
            count++;
            continue;
        }
        else
        {
            symbol = strtok (text, " =\n\r");

            assign_atom_parameters("symbol", symbol, (*atoms + count -2));

            for (i=0; i < DIM; i++)
            {
                char *convert;
                ccoord = strtok (NULL, " =\n\r");
                convert = &ccoord[0];
                (*atoms + count-2)->coords[i] = strtod(ccoord, &convert)/BOHR;
                (*atoms + count-2)->charge = ZERO;
            }
            count++;
        }
    }

    fclose(data);
    
    return natoms; 
}

void write_cubefile_offset(char * filename, cube_t * cube, int offset)
{
    FILE * datei;
    int i, j;
    int k, l, m;

    char spacer[]="  ";

    if((datei = fopen(filename, "w")))
    {
        fprintf(datei, "%s", cube->title);
        fprintf(datei, "%s", cube->comment);

        fprintf(datei, "%s%3i", spacer, cube->natoms);

        real origin[DIM];
        for (i=0; i<DIM; i++)
        {
            origin[i] = cube->origin[i];
            for (j=0; j<DIM; j++)
            {
                origin[i] += cube->boxv[i][j]*offset;
            }
        }

        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%10.6f", spacer, origin[i]);
        }

        fprintf(datei, "\n");

        for (i=0; i<DIM; i++)
        {
            fprintf(datei, "%s%3i", spacer, cube->n[i]-2*offset);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->boxv[i][j]);
            }
            fprintf(datei, "\n");
        }

        for (i=0; i<cube->natoms; i++)
        {
            fprintf(datei, "%s%3i%s%10.6f", spacer, cube->atoms[i].number, spacer, cube->atoms[i].charge);
            for (j=0; j<DIM; j++)
            {
                fprintf(datei, "%s%10.6f", spacer, cube->atoms[i].coords[j]);
            }
            fprintf(datei, "\n");
        }

        j = 0;
        
        int n0, n1, n2;
        n0 = cube->n[0] - offset;
        n1 = cube->n[1] - offset;
        n2 = cube->n[2] - offset;
        int index;

        for ( k=offset; k<n0; k++)
            for ( l=offset; l<n1; l++)
                for ( m=offset; m<n2; m++)
                {
                    index = m + cube->n[2] * ( l + cube->n[1] * k);
                    fprintf(datei, "%s%12.5e", spacer, cube->voxels[index].data);

                    if (j % 6 == 5)
                    {
                        fprintf(datei, "\n");
                    }

                    j++;
                }

        fprintf(datei, "\n");
    }

    fclose(datei);
}

void allocate_cube(voxel_t **voxels, int nvoxels)
{
    *voxels = (voxel_t*) malloc(nvoxels * (sizeof(voxel_t)));
}

void allocate_atoms(atom_t **atoms, int natoms)
{
    *atoms = (atom_t*) malloc(natoms * (sizeof(atom_t)));
}

void parse_cmdline(input_t * inppar, char ** argv, int argc)
{
    char * dummy;

    static const char *optString = "i:c:t:v:r:b:p:s:o:l:m:e:f:x:k:d:u:gnwh?";// Il:o:vh?";

    static const struct option longOpts[] = {
        { "input", required_argument, NULL, 'i' },
        { "cubefile", required_argument, NULL, 'c' },
        { "task", required_argument, NULL, 't' },
        { "verbose", required_argument, NULL, 'v' },
        { "reference", required_argument, NULL, 'r' },
        { "batch", required_argument, NULL, 'b' },
        { "prefix", required_argument, NULL, 'p' },
        { "suffix", required_argument, NULL, 's' },
        { "oprefix", required_argument, NULL, 'o' },
        { "lowbrow", required_argument, NULL, 'l' },
        { "mask", required_argument, NULL, 'm' },
        { "refmask", required_argument, NULL, 'e' },
        { "fudge", required_argument, NULL, 'f' },
        { "xmol", required_argument, NULL, 'x' },
        { "kurve", required_argument, NULL, 'k' },
        { "duplicate", required_argument, NULL, 'd' },
        { "resolution", required_argument, NULL, 'u' },
        { "guessfragments", no_argument, NULL, 'g' },
        { "nocore", no_argument, NULL, 'n' },
        { "wrap", no_argument, NULL, 'w' },
        // { "verbose", no_argument, NULL, 'v' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
    };

    int opt = 0;
    int longIndex = 0;

    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );

    while( opt != -1 ) 
    {
        switch( opt )
        {
            case 'i':
                strcpy(inppar->inpfile, optarg);
                break;
    
            case 't':
                strcpy(inppar->task, optarg);
                inppar-> tasknum = assign_task ( inppar->task );
                break;
    
            case 'v':
                inppar->output = atoi(optarg);
                break;
    
            case 'b':
                inppar->batchmode = 1;
    
                dummy = strtok(optarg, ":");
                inppar->start = atoi(dummy);
    
                dummy = strtok(NULL, ":");
                inppar->stop = atoi(dummy);
    
                dummy = strtok(NULL, ":");
                inppar->stride = atoi(dummy);
    
                break;
    
            case 'o':
                strcpy(inppar->outputprefix, optarg);
                break;
    
            case 'd':
                inppar->periodic = 1;
                /* check here if that is indeed okay, when we are working with cube files, but i see no problem so far */
                inppar->pbcset = 1;
    
                real conv = 1.;
                dummy = strtok(optarg, " ");

                if ( strstr ( dummy, "ANG" ) != NULL ) {
                    conv = ANG2BOHR;
                    dummy = strtok ( NULL, ":" );
                }
                else if ( strstr ( dummy, "NM" ) != NULL ) {
                    conv = NM2BOHR;
                    dummy = strtok ( NULL, ":" );
                }

                if ( dummy != NULL )
                    inppar->pbc[0] = conv * atof(dummy);
                else {
                        print_error ( MISSING_INPUT_PARAM, "cell length(s).");
                        exit ( MISSING_INPUT_PARAM );
                }
    
                dummy = strtok(NULL, ":");

                if ( dummy != NULL ) {
                    inppar->pbc[1] = conv * atof(dummy);

                    dummy = strtok(NULL, ":");
                    if ( dummy != NULL )
                        inppar->pbc[2] = conv * atof(dummy);
                    else {
                        print_error ( MISSING_INPUT_PARAM, "third periodic dimension. You only provided two dimensions for a 3D system.");
                        exit ( MISSING_INPUT_PARAM );
                    }

                }
                else {
                    inppar->pbc[1] = inppar->pbc[0];
                    inppar->pbc[2] = inppar->pbc[0];

                    break;
                }

                break;
            case 'u':
                dummy = strtok ( optarg, " " );

                inppar->resolution = atof ( dummy );
    
                break;
            case 'g':
                inppar->guessfragments = 1;
                break;
            case 'w':
                inppar->wrap = 1;
                break;
            case 'm':
                dummy = strtok(optarg, " ");
                strcpy(inppar->maskkind, dummy);

                dummy = strtok(NULL, " ");
                inppar->nkinds = 0;

                strcpy(inppar->mask, dummy);

                while ( dummy != NULL )
                {
                    dummy = strtok(NULL, " ");
                    if ( dummy != NULL)
                    {
                        strcat(inppar->mask, ",");
                        strcat(inppar->mask, dummy);
                    }
                    inppar->nkinds++;
                }
                break;
            case 'e':
                dummy = strtok(optarg, " ");
                strcpy(inppar->refmaskkind, dummy);

                dummy = strtok(NULL, " ");
                inppar->refnkinds = 0;

                strcpy(inppar->refmask, dummy);

                while ( dummy != NULL )
                {
                    dummy = strtok(NULL, " ");
                    if ( dummy != NULL)
                    {
                        strcat(inppar->refmask, ",");
                        strcat(inppar->refmask, dummy);
                    }
                    inppar->refnkinds++;
                }
                break;
            case 'f':
                inppar->blfudge = atof(optarg);
                break;
            case 'x':
                strcpy(inppar->structure, optarg);
                break;
            case 'k':
                strcpy(inppar->trajectory, optarg);

                if ( strstr(inppar->trajectory, ".xtc") != NULL )
#ifdef XDRCAP
                    inppar->xdrread = 1;
#else
                {
                    printf("SURF not compiled with XDR support\n");
                    exit ( 1 );
                }
#endif

                inppar->trajmode = 1;
                break;

            case 'h':   /* fall-through is intentional */
            case '?':
                display_usage();
                break;
    
            default:
                /* You won't actually get here. */
                // display_usage();
            break;
        }

        opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    }

    if ( optind != argc )
    {
        strcpy(inppar->structure, optarg);
    }
    
    return;
}

void display_usage(void)
{

    // printf("This is how it should be done...\n");
    printf("\nUsage: ccube [options] file...\n\n");
    printf("OPTIONS:\n");
    printf("    -i, --input=INPFILE             provide input file with all relevant data, options also given in command line will be overwritten\n");
    printf("    -t, --task=TASK                 task to be perfomed on cubefile(s)\n");
    printf("    -v, --verbose[=INTEGER]         output verbosity (0 - 3) \n");
    printf("    -b, --batch=START:STOP:STRIDE   enter batch processing mode; provide 3 integers separated by colons to determine start, stop and stride;\n");
    printf("                                    required additional input are --prefix and --suffix\n");
    printf("    -o, --oprefix=OPREFIX           prefix for output files\n");
    printf("    -m, --mask=MASK                 provide input mask (in apostrophes)\n");
    printf("    -e, --refmask=REFMASK           provide reference mask (in apostrophes)\n");
    printf("    -f, --fudge=FUDGE               fudge factor for acceptable bond lengths, i.e., multiplicative factor for sum of covalent radii\n");
    printf("    -x, --xmol=XMOL                 structure in xmol format, if provided overrides structure given in cube file\n");
    printf("    -k, --kurve=KURVE               trajectory file in XMOL format\n");
    printf("    -d, --duplicate=PBC             treat structure as it were duplicated periodically in space (i.e. pbc), 1-3 args; 1 arg: 3D with same PBC, 2-3: 2-3D\n");
    printf("    -u, --resolution                provide resolution for radial and angular investigation (either only rad or both)\n");
    printf("    -g, --guessfragments            guess fragments for analysis (currently only for AIM analysis)\n");
    printf("    -w, --wrap                      wrap molecules into box around (0,0,0) before doing any analysis\n");
    printf("    -h, --help                      show this message\n");
    printf("\n");
    
    exit(EXIT_FAILURE);
}

void xmolreader(FILE * fxmol, int bytelen, int snap, atom_t * atoms, int natoms)
{
    char text[MAXSTRLEN];
    char * dummy;
    int i, j;

    fseek ( fxmol, snap * bytelen, SEEK_SET );

    // get header information (which we HOPEFULLY) already obtained before
    fgets ( text, MAXSTRLEN, fxmol );
    fgets ( text, MAXSTRLEN, fxmol );

    for ( i=0; i<natoms; i++ )
    {
        fgets ( text, MAXSTRLEN, fxmol );

        dummy = strtok(text, " \n");

        for ( j=0; j<DIM; j++ )
        {
            dummy = strtok ( NULL, " \n");
            atoms[i].coords[j] = atof(dummy) / BOHR;
        }
    }

}

int xmol_snap_bytesize(FILE * fxmol)
{
    int natoms;
    int bytelen = 0;
    char text[MAXSTRLEN];

    rewind(fxmol);

    fgets ( text, MAXSTRLEN, fxmol );
    bytelen += strlen(text);
    natoms = atoi(text);

    fgets ( text, MAXSTRLEN, fxmol );
    bytelen += strlen(text);

    fgets ( text, MAXSTRLEN, fxmol );
    bytelen += ( natoms * strlen(text) );

    return bytelen; 
}

void read_xtr_forward ( XDRFILE * xd_read, int frwrd, atom_t * atoms, int natoms )
{
    int i, j;
    // XDRFILE *xd_read;
    rvec *x_xtc;
    float prec_xtc = 1000.0;
    float time_xtc;
    matrix box_xtc;
    int step_xtc;
    int result_xtc;
    
    // xd_read = xdrfile_open(rfile, "r");
    // result_xtc = read_xtc_natoms(rfile, &natoms_xtc);

    // if (exdrOK != result_xtc) {
    //     printf("Something went wrong opening XDR file\n"); // Error
    //     exit ( 1 );
    // }

    x_xtc = calloc(natoms, sizeof(x_xtc[0]));

    for ( i=0; i<frwrd; i++ )
#ifdef XDRCAP
        result_xtc = read_xtc(xd_read, natoms, &step_xtc, &time_xtc, box_xtc, x_xtc, &prec_xtc);
#endif

    if ( exdrOK != result_xtc ) {
        printf("Something went wrong reading the XTC trajectory.\nMaybe not enough snapshots?\n");
        exit ( 1 );
    }

    for ( i=0; i<natoms; i++ ) {
        // printf("%10.4f %10.4f %10.4f\n", x_xtc[i][0], x_xtc[i][1], x_xtc[i][2]);
        // printf("%i\n", i);
        for ( j=0; j<DIM; j++ )
            atoms[i].coords[j] = x_xtc[i][j] * NM2BOHR;
    }

    free ( x_xtc );
    // convert to useable format (for our purpose
}

void write_xyz ( FILE * xyz, atom_t * atoms, int natoms, char * comment )
{
    int i;

    fprintf(xyz, "%8i\n%s\n", natoms, comment);

    for ( i=0; i<natoms; i++ )
        fprintf(xyz, "%3s%21.10f%21.10f%21.10f\n", atoms[i].symbol, atoms[i].coords[0]*BOHR, atoms[i].coords[1]*BOHR, atoms[i].coords[2]*BOHR);
}

void write_matrix_real_2d_to_file_w_cont_spacing ( char * fname, real ** mat, int n1, int n2, real *origin, real d1, real d2 )
{
    FILE * datei;
    int i, j;

    if ( ( datei = fopen(fname, "w") ) )
    {
        fprintf(datei, "%14.8f", UNOBTAINABLE);

        for ( j=0; j<n2; j++ )
            fprintf ( datei, " %14.8f", origin[1]+d2*j );

        fprintf(datei, "\n");

        for ( i=0; i<n1; i++ ) {
            fprintf ( datei, " %14.8f", origin[0] + d1*i);
            for ( j=0; j<n2; j++ )
                fprintf(datei, " %14.8f", mat[i][j]);

            fprintf(datei, "\n");
        }

        fclose ( datei );
    }
    else {
        printf("Could not open roots file '%s'\n", fname);
        exit ( 1 );
    }
}

int assign_task ( char * task )
{
    if ( strstr ( task, "surface_distribution" ) != NULL )
        return SURFDIST;
    else if ( strstr ( task, "topology" ) != NULL )
        return TOPOLOGY;
    else
        return SURFDIST;
}
