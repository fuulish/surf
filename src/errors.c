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
#include "constants.h"
#include "atom_param.h"
#include <stdio.h>
#include <string.h>

void print_error(int errno, char * detail)
{

    print_error_header();

    switch ( errno )
    {
        case FILE_NOT_FOUND:
            printf("File '%s' could not be opened.\n", detail);
            break;
        case MISSING_INPUT_PARAM:
            printf("Task cannot be completed without specifying input parameter: %s\n", detail);
            break;
        case NOT_IMPLEMENTED:
            printf("%s is not implemented, yet.\n", detail);
            break;
        case PROGRAM_BROKEN:
            printf("Trying %s will break the program. A workaround may be implemented soon. Stay tuned!\n", detail);
            break;
        case OUT_OF_MEMORY:
            printf("Out of memory during %s.\n", detail);
            break;
        default:
            printf("An error occurred that has not yet been assigned to a category!\n");
            break;
    }

    print_error_footer();
}

void print_error_header()
{
    printf("\n!!!!!!!!! E R R O R !!!!!!!!!\n\n");
}

void print_error_footer()
{
    printf("\n!!!!!!!!! E R R O R !!!!!!!!!\n\n");
}
