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
#include "errors.h"
#include "io.h"
#include "constants.h"
#include "cube.h"
#include "molmanipul.h"
#include "trajanal.h"
#include "main.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main (int argc, char *argv[])
{
    input_t inppar;
    int t;

    // printf("You are running SURF version SURF_VERSION_MAJOR.SURF_VERSION_MINOR\n");
    set_input_defaults(&inppar);

    if (argc == 1)
    {
        display_usage();
    }
    else if (argc == 2)
    {
        if ( strncmp ( argv[1], "-", 1 )  == 0) // (strncmp ( argv[1], "-h", 2 ) == 0 ) || ( strncmp ( argv[1], "--help", 5) == 0 ) )
        {
            display_usage();
        }
        else
        {
            strcpy(inppar.structure, argv[1]);
        }
    }
    else
    {
        parse_cmdline(&inppar, argv, argc);
    }

    if (strncmp(inppar.inpfile, EMPTY, 10) != 0)
        parse_input_file(&inppar, inppar.inpfile);

    if  ( ( inppar.batchmode ) && ( inppar.trajmode ) ) {
        tanalize(&inppar);
    }

    t = clock();

    printf("-----------------------------------------\n\n");
    printf("Total cpu time:      %4.2f s\n", ((float)t)/CLOCKS_PER_SEC);
    printf("-----------------------------------------\n\n");

    return 0;
}

