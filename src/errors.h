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

typedef enum errors_e
{
    FILE_NOT_FOUND=100,
    MISSING_INPUT_PARAM,
    PROGRAM_BROKEN,
    OUT_OF_MEMORY,
    MISSING_LIBRARY,
    CONFLICTING_OPTIONS,
    NOT_IMPLEMENTED,

} errors_t;

void print_error(int errno, char * detail);
void print_error_header();
void print_error_footer();
