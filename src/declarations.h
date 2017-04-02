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

#define DIM 3
#define sqr(x) ((x)*(x))
#define FLEN sizeof(float)
#define CLEN sizeof(char[1])
#define DLLEN sizeof(long real)
#define NUMKEYS 36
#define EPS 1.e-6
#define ZERO 0.
#define INTZERO 0
#define ONE 1.
#define EMPTY "empty like a stallion's pussy"
#define MAXSTRLEN 100000
#define UNAVAIL 0.
#define UNOBTAINABLE -2468.13579
#define HELLO printf("HELLO THERE!\n")
