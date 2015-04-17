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
#include <stdlib.h>

void assign_atom_parameters( char * kind, char * symdex , atom_t * atom)
{
    // int index;
    // otherwise it maybe unitialized and the code would break, check here again later
    int index = ATOM_DUMM;

    if (strncmp(kind, "index", 5) == 0)
    {
        index = atoi(symdex);
    }
    if (strncmp(kind, "symbol", 5) == 0)
    {
        if (strncmp(symdex, "H", 2) == 0)
        {
            index = ATOM_H;
        }
        else if (strncmp(symdex, "He", 2) == 0)
        {
            index = ATOM_HE;
        }
        else if (strncmp(symdex, "Li", 2) == 0)
        {
            index = ATOM_LI;
        }
        else if (strncmp(symdex, "Be", 2) == 0)
        {
            index = ATOM_BE;
        }
        else if (strncmp(symdex, "B", 2) == 0)
        {
            index = ATOM_B;
        }
        else if (strncmp(symdex, "C", 2) == 0)
        {
            index = ATOM_C;
        }
        else if (strncmp(symdex, "N", 2) == 0)
        {
            index = ATOM_N;
        }
        else if (strncmp(symdex, "O", 2) == 0)
        {
            index = ATOM_O;
        }
        else if (strncmp(symdex, "F", 2) == 0)
        {
            index = ATOM_F;
        }
        else if (strncmp(symdex, "Ne", 2) == 0)
        {
            index = ATOM_NE;
        }
        else if (strncmp(symdex, "Na", 2) == 0)
        {
            index = ATOM_NA;
        }
        else if (strncmp(symdex, "Mg", 2) == 0)
        {
            index = ATOM_MG;
        }
        else if (strncmp(symdex, "Al", 2) == 0)
        {
            index = ATOM_AL;
        }
        else if (strncmp(symdex, "Si", 2) == 0)
        {
            index = ATOM_SI;
        }
        else if (strncmp(symdex, "P", 2) == 0)
        {
            index = ATOM_P;
        }
        else if (strncmp(symdex, "S", 2) == 0)
        {
            index = ATOM_S;
        }
        else if (strncmp(symdex, "Cl", 2) == 0)
        {
            index = ATOM_CL;
        }
        else if (strncmp(symdex, "Ar", 2) == 0)
        {
            index = ATOM_AR;
        }
        else if (strncmp(symdex, "K", 2) == 0)
        {
            index = ATOM_K;
        }
        else if (strncmp(symdex, "Ca", 2) == 0)
        {
            index = ATOM_CA;
        }
        else if (strncmp(symdex, "Sc", 2) == 0)
        {
            index = ATOM_SC;
        }
        else if (strncmp(symdex, "Ti", 2) == 0)
        {
            index = ATOM_TI;
        }
        else if (strncmp(symdex, "V", 2) == 0)
        {
            index = ATOM_V;
        }
        else if (strncmp(symdex, "Cr", 2) == 0)
        {
            index = ATOM_CR;
        }
        else if (strncmp(symdex, "Mn", 2) == 0)
        {
            index = ATOM_MN;
        }
        else if (strncmp(symdex, "Fe", 2) == 0)
        {
            index = ATOM_FE;
        }
        else if (strncmp(symdex, "Co", 2) == 0)
        {
            index = ATOM_CO;
        }
        else if (strncmp(symdex, "Ni", 2) == 0)
        {
            index = ATOM_NI;
        }
        else if (strncmp(symdex, "Cu", 2) == 0)
        {
            index = ATOM_CU;
        }
        else if (strncmp(symdex, "Zn", 2) == 0)
        {
            index = ATOM_ZN;
        }
        else if (strncmp(symdex, "Ga", 2) == 0)
        {
            index = ATOM_GA;
        }
        else if (strncmp(symdex, "Ge", 2) == 0)
        {
            index = ATOM_GE;
        }
        else if (strncmp(symdex, "As", 2) == 0)
        {
            index = ATOM_AS;
        }
        else if (strncmp(symdex, "Se", 2) == 0)
        {
            index = ATOM_SE;
        }
        else if (strncmp(symdex, "Br", 2) == 0)
        {
            index = ATOM_BR;
        }
        else if (strncmp(symdex, "Kr", 2) == 0)
        {
            index = ATOM_KR;
        }
        else if (strncmp(symdex, "Rb", 2) == 0)
        {
            index = ATOM_RB;
        }
        else if (strncmp(symdex, "Sr", 2) == 0)
        {
            index = ATOM_SR;
        }
        else if (strncmp(symdex, "Y", 2) == 0)
        {
            index = ATOM_Y;
        }
        else if (strncmp(symdex, "Zr", 2) == 0)
        {
            index = ATOM_ZR;
        }
        else if (strncmp(symdex, "Nb", 2) == 0)
        {
            index = ATOM_NB;
        }
        else if (strncmp(symdex, "Mo", 2) == 0)
        {
            index = ATOM_MO;
        }
        else if (strncmp(symdex, "Tc", 2) == 0)
        {
            index = ATOM_TC;
        }
        else if (strncmp(symdex, "Ru", 2) == 0)
        {
            index = ATOM_RU;
        }
        else if (strncmp(symdex, "Rh", 2) == 0)
        {
            index = ATOM_RH;
        }
        else if (strncmp(symdex, "Pd", 2) == 0)
        {
            index = ATOM_PD;
        }
        else if (strncmp(symdex, "Ag", 2) == 0)
        {
            index = ATOM_AG;
        }
        else if (strncmp(symdex, "Cd", 2) == 0)
        {
            index = ATOM_CD;
        }
        else if (strncmp(symdex, "In", 2) == 0)
        {
            index = ATOM_IN;
        }
        else if (strncmp(symdex, "Sn", 2) == 0)
        {
            index = ATOM_SN;
        }
        else if (strncmp(symdex, "Sb", 2) == 0)
        {
            index = ATOM_SB;
        }
        else if (strncmp(symdex, "Te", 2) == 0)
        {
            index = ATOM_TE;
        }
        else if (strncmp(symdex, "I", 2) == 0)
        {
            index = ATOM_I;
        }
        else if (strncmp(symdex, "Xe", 2) == 0)
        {
            index = ATOM_XE;
        }
        else if (strncmp(symdex, "Cs", 2) == 0)
        {
            index = ATOM_CS;
        }
        else if (strncmp(symdex, "Ba", 2) == 0)
        {
            index = ATOM_BA;
        }
        else if (strncmp(symdex, "La", 2) == 0)
        {
            index = ATOM_LA;
        }
        else if (strncmp(symdex, "Ce", 2) == 0)
        {
            index = ATOM_CE;
        }
        else if (strncmp(symdex, "Pr", 2) == 0)
        {
            index = ATOM_PR;
        }
        else if (strncmp(symdex, "Nd", 2) == 0)
        {
            index = ATOM_ND;
        }
        else if (strncmp(symdex, "Pm", 2) == 0)
        {
            index = ATOM_PM;
        }
        else if (strncmp(symdex, "Sm", 2) == 0)
        {
            index = ATOM_SM;
        }
        else if (strncmp(symdex, "Eu", 2) == 0)
        {
            index = ATOM_EU;
        }
        else if (strncmp(symdex, "Gd", 2) == 0)
        {
            index = ATOM_GD;
        }
        else if (strncmp(symdex, "Tb", 2) == 0)
        {
            index = ATOM_TB;
        }
        else if (strncmp(symdex, "Dy", 2) == 0)
        {
            index = ATOM_DY;
        }
        else if (strncmp(symdex, "Ho", 2) == 0)
        {
            index = ATOM_HO;
        }
        else if (strncmp(symdex, "Er", 2) == 0)
        {
            index = ATOM_ER;
        }
        else if (strncmp(symdex, "Tm", 2) == 0)
        {
            index = ATOM_TM;
        }
        else if (strncmp(symdex, "Yb", 2) == 0)
        {
            index = ATOM_YB;
        }
        else if (strncmp(symdex, "Lu", 2) == 0)
        {
            index = ATOM_LU;
        }
        else if (strncmp(symdex, "Hf", 2) == 0)
        {
            index = ATOM_HF;
        }
        else if (strncmp(symdex, "Ta", 2) == 0)
        {
            index = ATOM_TA;
        }
        else if (strncmp(symdex, "W", 2) == 0)
        {
            index = ATOM_W;
        }
        else if (strncmp(symdex, "Re", 2) == 0)
        {
            index = ATOM_RE;
        }
        else if (strncmp(symdex, "Os", 2) == 0)
        {
            index = ATOM_OS;
        }
        else if (strncmp(symdex, "Ir", 2) == 0)
        {
            index = ATOM_IR;
        }
        else if (strncmp(symdex, "Pt", 2) == 0)
        {
            index = ATOM_PT;
        }
        else if (strncmp(symdex, "Au", 2) == 0)
        {
            index = ATOM_AU;
        }
        else if (strncmp(symdex, "Hg", 2) == 0)
        {
            index = ATOM_HG;
        }
        else if (strncmp(symdex, "Tl", 2) == 0)
        {
            index = ATOM_TL;
        }
        else if (strncmp(symdex, "Pb", 2) == 0)
        {
            index = ATOM_PB;
        }
        else if (strncmp(symdex, "Bi", 2) == 0)
        {
            index = ATOM_BI;
        }
        else if (strncmp(symdex, "Po", 2) == 0)
        {
            index = ATOM_PO;
        }
        else if (strncmp(symdex, "At", 2) == 0)
        {
            index = ATOM_AT;
        }
        else if (strncmp(symdex, "Rn", 2) == 0)
        {
            index = ATOM_RN;
        }
        else if (strncmp(symdex, "Fr", 2) == 0)
        {
            index = ATOM_FR;
        }
        else if (strncmp(symdex, "Ra", 2) == 0)
        {
            index = ATOM_RA;
        }
        else if (strncmp(symdex, "Ac", 2) == 0)
        {
            index = ATOM_AC;
        }
        else if (strncmp(symdex, "Th", 2) == 0)
        {
            index = ATOM_TH;
        }
        else if (strncmp(symdex, "Pa", 2) == 0)
        {
            index = ATOM_PA;
        }
        else if (strncmp(symdex, "U", 2) == 0)
        {
            index = ATOM_U;
        }
        else if (strncmp(symdex, "Np", 2) == 0)
        {
            index = ATOM_NP;
        }
        else if (strncmp(symdex, "Pu", 2) == 0)
        {
            index = ATOM_PU;
        }
        else if (strncmp(symdex, "Am", 2) == 0)
        {
            index = ATOM_AM;
        }
        else if (strncmp(symdex, "Cm", 2) == 0)
        {
            index = ATOM_CM;
        }
        else if (strncmp(symdex, "Bk", 2) == 0)
        {
            index = ATOM_BK;
        }
        else if (strncmp(symdex, "Cf", 2) == 0)
        {
            index = ATOM_CF;
        }
        else if (strncmp(symdex, "Es", 2) == 0)
        {
            index = ATOM_ES;
        }
        else if (strncmp(symdex, "Fm", 2) == 0)
        {
            index = ATOM_FM;
        }
        else if (strncmp(symdex, "Md", 2) == 0)
        {
            index = ATOM_MD;
        }
        else if (strncmp(symdex, "No", 2) == 0)
        {
            index = ATOM_NO;
        }
        else if (strncmp(symdex, "Lr", 2) == 0)
        {
            index = ATOM_LR;
        }
        else if (strncmp(symdex, "Rf", 2) == 0)
        {
            index = ATOM_RF;
        }
        else if (strncmp(symdex, "Db", 2) == 0)
        {
            index = ATOM_DB;
        }
        else if (strncmp(symdex, "Sg", 2) == 0)
        {
            index = ATOM_SG;
        }
        else if (strncmp(symdex, "Bh", 2) == 0)
        {
            index = ATOM_BH;
        }
        else if (strncmp(symdex, "Hs", 2) == 0)
        {
            index = ATOM_HS;
        }
        else if (strncmp(symdex, "Mt", 2) == 0)
        {
            index = ATOM_MT;
        }
        else if (strncmp(symdex, "Ds", 2) == 0)
        {
            index = ATOM_DS;
        }
        else if (strncmp(symdex, "Rg", 2) == 0)
        {
            index = ATOM_RG;
        }
        else if (strncmp(symdex, "Cn", 2) == 0)
        {
            index = ATOM_CN;
        }
        else if (strncmp(symdex, "Uut", 2) == 0)
        {
            index = ATOM_UUT;
        }
        else if (strncmp(symdex, "Fl", 2) == 0)
        {
            index = ATOM_FL;
        }
        else if (strncmp(symdex, "Uup", 2) == 0)
        {
            index = ATOM_UUP;
        }
        else if (strncmp(symdex, "Lv", 2) == 0)
        {
            index = ATOM_LV;
        }
        else if (strncmp(symdex, "Uus", 2) == 0)
        {
            index = ATOM_UUS;
        }
        else if (strncmp(symdex, "Uuo", 2) == 0)
        {
            index = ATOM_UUO;
        }
        else if (strncmp(symdex, "Du", 2) == 0)
        {
            index = ATOM_DUMM;
        }
        else if (strncmp(symdex, "X", 2) == 0)
        {
            index = ATOM_WANN;
        }
    }

    switch ( index )
    {
        case ATOM_H:
            atom->covrad = 37/100./BOHR;
            atom->rvdw = 120./100./BOHR;
            atom->charge = ATOM_H;
            atom->number = ATOM_H;
            strcpy(atom->symbol, "H");
            atom->coreel = 0;
            atom->mass = 1.0079;
            break;
        case ATOM_HE:
            atom->covrad = 32/100./BOHR;
            atom->rvdw = 140./100./BOHR;
            atom->charge = ATOM_HE;
            atom->number = ATOM_HE;
            strcpy(atom->symbol, "He");
            atom->coreel = 0;
            atom->mass = 4.0026;
            break;
        case ATOM_LI:
            atom->covrad = 134/100./BOHR;
            atom->rvdw = 182./100./BOHR;
            atom->charge = ATOM_LI;
            atom->number = ATOM_LI;
            strcpy(atom->symbol, "Li");
            atom->coreel = 2;
            atom->mass = 6.941;
            break;
        case ATOM_BE:
            atom->covrad = 90/100./BOHR;
            atom->rvdw = 153./100./BOHR;
            atom->charge = ATOM_BE;
            atom->number = ATOM_BE;
            strcpy(atom->symbol, "Be");
            atom->coreel = 2;
            atom->mass = 9.0122;
            break;
        case ATOM_B:
            atom->covrad = 82/100./BOHR;
            atom->rvdw = 192./100./BOHR;
            atom->charge = ATOM_B;
            atom->number = ATOM_B;
            strcpy(atom->symbol, "B");
            atom->coreel = 2;
            atom->mass = 10.811;
            break;
        case ATOM_C:
            atom->covrad = 77/100./BOHR;
            atom->rvdw = 170./100./BOHR;
            atom->charge = ATOM_C;
            atom->number = ATOM_C;
            strcpy(atom->symbol, "C");
            atom->coreel = 2;
            atom->mass = 12.0107;
            break;
        case ATOM_N:
            atom->covrad = 75/100./BOHR;
            atom->rvdw = 155./100./BOHR;
            atom->charge = ATOM_N;
            atom->number = ATOM_N;
            strcpy(atom->symbol, "N");
            atom->coreel = 2;
            atom->mass = 14.0067;
            break;
        case ATOM_O:
            atom->covrad = 73/100./BOHR;
            atom->rvdw = 152./100./BOHR;
            atom->charge = ATOM_O;
            atom->number = ATOM_O;
            strcpy(atom->symbol, "O");
            atom->coreel = 2;
            atom->mass = 15.9994;
            break;
        case ATOM_F:
            atom->covrad = 71/100./BOHR;
            atom->rvdw = 147./100./BOHR;
            atom->charge = ATOM_F;
            atom->number = ATOM_F;
            strcpy(atom->symbol, "F");
            atom->coreel = 2;
            atom->mass = 18.9984;
            break;
        case ATOM_NE:
            atom->covrad = 69/100./BOHR;
            atom->rvdw = 154./100./BOHR;
            atom->charge = ATOM_NE;
            atom->number = ATOM_NE;
            strcpy(atom->symbol, "Ne");
            atom->coreel = 2;
            atom->mass = 20.1797;
            break;
        case ATOM_NA:
            atom->covrad = 154/100./BOHR;
            atom->rvdw = 227./100./BOHR;
            atom->charge = ATOM_NA;
            atom->number = ATOM_NA;
            strcpy(atom->symbol, "Na");
            atom->coreel = 2;
            atom->mass = 22.9897;
            break;
        case ATOM_MG:
            atom->covrad = 130/100./BOHR;
            atom->rvdw = 173./100./BOHR;
            atom->charge = ATOM_MG;
            atom->number = ATOM_MG;
            strcpy(atom->symbol, "Mg");
            atom->coreel = 10;
            atom->mass = 24.305;
            break;
        case ATOM_AL:
            atom->covrad = 118/100./BOHR;
            atom->rvdw = 184./100./BOHR;
            atom->charge = ATOM_AL;
            atom->number = ATOM_AL;
            strcpy(atom->symbol, "Al");
            atom->coreel = 10;
            atom->mass = 26.9815;
            break;
        case ATOM_SI:
            atom->covrad = 111/100./BOHR;
            atom->rvdw = 210./100./BOHR;
            atom->charge = ATOM_SI;
            atom->number = ATOM_SI;
            strcpy(atom->symbol, "Si");
            atom->coreel = 10;
            atom->mass = 28.0855;
            break;
        case ATOM_P:
            atom->covrad = 106/100./BOHR;
            atom->rvdw = 180./100./BOHR;
            atom->charge = ATOM_P;
            atom->number = ATOM_P;
            strcpy(atom->symbol, "P");
            atom->coreel = 10;
            atom->mass = 30.9738;
            break;
        case ATOM_S:
            atom->covrad = 102/100./BOHR;
            atom->rvdw = 180./100./BOHR;
            atom->charge = ATOM_S;
            atom->number = ATOM_S;
            strcpy(atom->symbol, "S");
            atom->coreel = 10;
            atom->mass = 32.065;
            break;
        case ATOM_CL:
            atom->covrad = 99/100./BOHR;
            atom->rvdw = 175./100./BOHR;
            atom->charge = ATOM_CL;
            atom->number = ATOM_CL;
            strcpy(atom->symbol, "Cl");
            atom->coreel = 10;
            atom->mass = 35.453;
            break;
        case ATOM_AR:
            atom->covrad = 97/100./BOHR;
            atom->rvdw = 188./100./BOHR;
            atom->charge = ATOM_AR;
            atom->number = ATOM_AR;
            strcpy(atom->symbol, "Ar");
            atom->coreel = 10;
            atom->mass = 39.948;
            break;
        case ATOM_K:
            atom->covrad = 196/100./BOHR;
             atom->rvdw = 275./100./BOHR;
             atom->charge = ATOM_K;
             atom->number = ATOM_K;
             strcpy(atom->symbol, "K");
             atom->coreel = 18;
            atom->mass = 39.0983;
             break;
        case ATOM_CA:
            atom->covrad = 174/100./BOHR;
            atom->rvdw = 231./100./BOHR;
            atom->charge = ATOM_CA;
            atom->number = ATOM_CA;
            strcpy(atom->symbol, "Ca");
            atom->coreel = 18;
            atom->mass = 40.078;
            break;
        case ATOM_SC:
            atom->covrad = 144/100./BOHR;
            atom->rvdw = 211./100./BOHR;
            atom->charge = ATOM_SC;
            atom->number = ATOM_SC;
            strcpy(atom->symbol, "Sc");
            atom->coreel = 18;
            atom->mass = 44.9559;
            break;
        case ATOM_TI:
            atom->covrad = 136/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_TI;
            atom->number = ATOM_TI;
            strcpy(atom->symbol, "Ti");
            atom->coreel = 18;
            atom->mass = 47.867;
            break;
        case ATOM_V:
            atom->covrad = 125/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_V;
            atom->number = ATOM_V;
            strcpy(atom->symbol, "V");
            atom->coreel = 18;
            atom->mass = 50.9415;
            break;
        case ATOM_CR:
            atom->covrad = 127/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_CR;
            atom->number = ATOM_CR;
            strcpy(atom->symbol, "Cr");
            atom->coreel = 18;
            atom->mass = 51.9961;
            break;
        case ATOM_MN:
            atom->covrad = 139/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_MN;
            atom->number = ATOM_MN;
            strcpy(atom->symbol, "Mn");
            atom->coreel = 18;
            atom->mass = 54.938;
            break;
        case ATOM_FE:
            atom->covrad = 125/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_FE;
            atom->number = ATOM_FE;
            strcpy(atom->symbol, "Fe");
            atom->coreel = 18;
            atom->mass = 55.845;
            break;
        case ATOM_CO:
            atom->covrad = 126/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_CO;
            atom->number = ATOM_CO;
            strcpy(atom->symbol, "Co");
            atom->coreel = 18;
            atom->mass = 58.9332;
            break;
        case ATOM_NI:
            atom->covrad = 121/100./BOHR;
            atom->rvdw = 163./100./BOHR;
            atom->charge = ATOM_NI;
            atom->number = ATOM_NI;
            strcpy(atom->symbol, "Ni");
            atom->coreel = 18;
            atom->mass = 58.6934;
            break;
        case ATOM_CU:
            atom->covrad = 138/100./BOHR;
            atom->rvdw = 140./100./BOHR;
            atom->charge = ATOM_CU;
            atom->number = ATOM_CU;
            strcpy(atom->symbol, "Cu");
            atom->coreel = 18;
            atom->mass = 63.546;
            break;
        case ATOM_ZN:
            atom->covrad = 131/100./BOHR;
            atom->rvdw = 139./100./BOHR;
            atom->charge = ATOM_ZN;
            atom->number = ATOM_ZN;
            strcpy(atom->symbol, "Zn");
            atom->coreel = 18;
            atom->mass = 65.39;
            break;
        case ATOM_GA:
            atom->covrad = 126/100./BOHR;
            atom->rvdw = 187./100./BOHR;
            atom->charge = ATOM_GA;
            atom->number = ATOM_GA;
            strcpy(atom->symbol, "Ga");
            atom->coreel = 18;
            atom->mass = 69.723;
            break;
        case ATOM_GE:
            atom->covrad = 122/100./BOHR;
            atom->rvdw = 211./100./BOHR;
            atom->charge = ATOM_GE;
            atom->number = ATOM_GE;
            strcpy(atom->symbol, "Ge");
            atom->coreel = 18;
            atom->mass = 72.64;
            break;
        case ATOM_AS:
            atom->covrad = 119/100./BOHR;
            atom->rvdw = 185./100./BOHR;
            atom->charge = ATOM_AS;
            atom->number = ATOM_AS;
            strcpy(atom->symbol, "As");
            atom->coreel = 18;
            atom->mass = 74.9216;
            break;
        case ATOM_SE:
            atom->covrad = 116/100./BOHR;
            atom->rvdw = 190./100./BOHR;
            atom->charge = ATOM_SE;
            atom->number = ATOM_SE;
            strcpy(atom->symbol, "Se");
            atom->coreel = 18;
            atom->mass = 78.96;
            break;
        case ATOM_BR:
            atom->covrad = 114/100./BOHR;
            atom->rvdw = 185./100./BOHR;
            atom->charge = ATOM_BR;
            atom->number = ATOM_BR;
            strcpy(atom->symbol, "Br");
            atom->coreel = 18;
            atom->mass = 79.904;
            break;
        case ATOM_KR:
            atom->covrad = 110/100./BOHR;
            atom->rvdw = 202./100./BOHR;
            atom->charge = ATOM_KR;
            atom->number = ATOM_KR;
            strcpy(atom->symbol, "Kr");
            atom->coreel = 18;
            atom->mass = 83.8;
            break;
        case ATOM_RB:
            atom->covrad = 211/100./BOHR;
            atom->rvdw = 303./100./BOHR;
            atom->charge = ATOM_RB;
            atom->number = ATOM_RB;
            strcpy(atom->symbol, "Rb");
            atom->coreel = 36;
            atom->mass = 85.4678;
            break;
        case ATOM_SR:
            atom->covrad = 192/100./BOHR;
            atom->rvdw = 249./100./BOHR;
            atom->charge = ATOM_SR;
            atom->number = ATOM_SR;
            strcpy(atom->symbol, "Sr");
            atom->coreel = 36;
            atom->mass = 87.62;
            break;
        case ATOM_Y:
            atom->covrad = 162/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_Y;
            atom->number = ATOM_Y;
            strcpy(atom->symbol, "Y");
            atom->coreel = 36;
            atom->mass = 88.9059;
            break;
        case ATOM_ZR:
            atom->covrad = 148/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_ZR;
            atom->number = ATOM_ZR;
            strcpy(atom->symbol, "Zr");
            atom->coreel = 36;
            atom->mass = 91.224;
            break;
        case ATOM_NB:
            atom->covrad = 137/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_NB;
            atom->number = ATOM_NB;
            strcpy(atom->symbol, "Nb");
            atom->coreel = 36;
            atom->mass = 92.9064;
            break;
        case ATOM_MO:
            atom->covrad = 145/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_MO;
            atom->number = ATOM_MO;
            strcpy(atom->symbol, "Mo");
            atom->coreel = 36;
            atom->mass = 95.94;
            break;
        case ATOM_TC:
            atom->covrad = 156/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_TC;
            atom->number = ATOM_TC;
            strcpy(atom->symbol, "Tc");
            atom->coreel = 36;
            atom->mass = 98;
            break;
        case ATOM_RU:
            atom->covrad = 126/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_RU;
            atom->number = ATOM_RU;
            strcpy(atom->symbol, "Ru");
            atom->coreel = 36;
            atom->mass = 101.07;
            break;
        case ATOM_RH:
            atom->covrad = 135/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_RH;
            atom->number = ATOM_RH;
            strcpy(atom->symbol, "Rh");
            atom->coreel = 36;
            atom->mass = 102.9055;
            break;
        case ATOM_PD:
            atom->covrad = 131/100./BOHR;
            atom->rvdw = 163./100./BOHR;
            atom->charge = ATOM_PD;
            atom->number = ATOM_PD;
            strcpy(atom->symbol, "Pd");
            atom->coreel = 36;
            atom->mass = 106.42;
            break;
        case ATOM_AG:
            atom->covrad = 153/100./BOHR;
            atom->rvdw = 172./100./BOHR;
            atom->charge = ATOM_AG;
            atom->number = ATOM_AG;
            strcpy(atom->symbol, "Ag");
            atom->coreel = 36;
            atom->mass = 107.8682;
            break;
        case ATOM_CD:
            atom->covrad = 148/100./BOHR;
            atom->rvdw = 158./100./BOHR;
            atom->charge = ATOM_CD;
            atom->number = ATOM_CD;
            strcpy(atom->symbol, "Cd");
            atom->coreel = 36;
            atom->mass = 112.411;
            break;
        case ATOM_IN:
            atom->covrad = 144/100./BOHR;
            atom->rvdw = 193./100./BOHR;
            atom->charge = ATOM_IN;
            atom->number = ATOM_IN;
            strcpy(atom->symbol, "In");
            atom->coreel = 36;
            atom->mass = 114.818;
            break;
        case ATOM_SN:
            atom->covrad = 141/100./BOHR;
            atom->rvdw = 217./100./BOHR;
            atom->charge = ATOM_SN;
            atom->number = ATOM_SN;
            strcpy(atom->symbol, "Sn");
            atom->coreel = 36;
            atom->mass = 118.71;
            break;
        case ATOM_SB:
            atom->covrad = 138/100./BOHR;
            atom->rvdw = 206./100./BOHR;
            atom->charge = ATOM_SB;
            atom->number = ATOM_SB;
            strcpy(atom->symbol, "Sb");
            atom->coreel = 36;
            atom->mass = 121.76;
            break;
        case ATOM_TE:
            atom->covrad = 135/100./BOHR;
            atom->rvdw = 206./100./BOHR;
            atom->charge = ATOM_TE;
            atom->number = ATOM_TE;
            strcpy(atom->symbol, "Te");
            atom->coreel = 36;
            atom->mass = 127.6;
            break;
        case ATOM_I:
            atom->covrad = 133/100./BOHR;
            atom->rvdw = 198./100./BOHR;
            atom->charge = ATOM_I;
            atom->number = ATOM_I;
            strcpy(atom->symbol, "I");
            atom->coreel = 36;
            atom->mass = 126.9045;
            break;
        case ATOM_XE:
            atom->covrad = 130/100./BOHR;
            atom->rvdw = 216./100./BOHR;
            atom->charge = ATOM_XE;
            atom->number = ATOM_XE;
            strcpy(atom->symbol, "Xe");
            atom->coreel = 36;
            atom->mass = 131.293;
            break;
        case ATOM_CS:
            atom->covrad = 225/100./BOHR;
            atom->rvdw = 343./100./BOHR;
            atom->charge = ATOM_CS;
            atom->number = ATOM_CS;
            strcpy(atom->symbol, "Cs");
            atom->coreel = 54;
            atom->mass = 132.9055;
            break;
        case ATOM_BA:
            atom->covrad = 198/100./BOHR;
            atom->rvdw = 268./100./BOHR;
            atom->charge = ATOM_BA;
            atom->number = ATOM_BA;
            strcpy(atom->symbol, "Ba");
            atom->coreel = 54;
            atom->mass = 137.327;
            break;
        case ATOM_LA:
            atom->covrad = 169/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_LA;
            atom->number = ATOM_LA;
            strcpy(atom->symbol, "La");
            atom->coreel = 54;
            atom->mass = 138.9055;
            break;
        case ATOM_CE:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_CE;
            atom->number = ATOM_CE;
            strcpy(atom->symbol, "Ce");
            atom->coreel = 54;
            atom->mass = 140.116;
            break;
        case ATOM_PR:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_PR;
            atom->number = ATOM_PR;
            strcpy(atom->symbol, "Pr");
            atom->coreel = 54;
            atom->mass = 140.9077;
            break;
        case ATOM_ND:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_ND;
            atom->number = ATOM_ND;
            strcpy(atom->symbol, "Nd");
            atom->coreel = 54;
            atom->mass = 144.24;
            break;
        case ATOM_PM:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_PM;
            atom->number = ATOM_PM;
            strcpy(atom->symbol, "Pm");
            atom->coreel = 54;
            atom->mass = 145;
            break;
        case ATOM_SM:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_SM;
            atom->number = ATOM_SM;
            strcpy(atom->symbol, "Sm");
            atom->coreel = 54;
            atom->mass = 150.36;
            break;
        case ATOM_EU:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_EU;
            atom->number = ATOM_EU;
            strcpy(atom->symbol, "Eu");
            atom->coreel = 54;
            atom->mass = 151.964;
            break;
        case ATOM_GD:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_GD;
            atom->number = ATOM_GD;
            strcpy(atom->symbol, "Gd");
            atom->coreel = 54;
            atom->mass = 157.25;
            break;
        case ATOM_TB:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_TB;
            atom->number = ATOM_TB;
            strcpy(atom->symbol, "Tb");
            atom->coreel = 54;
            atom->mass = 158.9253;
            break;
        case ATOM_DY:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_DY;
            atom->number = ATOM_DY;
            strcpy(atom->symbol, "Dy");
            atom->coreel = 54;
            atom->mass = 162.5;
            break;
        case ATOM_HO:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_HO;
            atom->number = ATOM_HO;
            strcpy(atom->symbol, "Ho");
            atom->coreel = 54;
            atom->mass = 164.9303;
            break;
        case ATOM_ER:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_ER;
            atom->number = ATOM_ER;
            strcpy(atom->symbol, "Er");
            atom->coreel = 54;
            atom->mass = 167.259;
            break;
        case ATOM_TM:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_TM;
            atom->number = ATOM_TM;
            strcpy(atom->symbol, "Tm");
            atom->coreel = 54;
            atom->mass = 168.9342;
            break;
        case ATOM_YB:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_YB;
            atom->number = ATOM_YB;
            strcpy(atom->symbol, "Yb");
            atom->coreel = 54;
            atom->mass = 173.04;
            break;
        case ATOM_LU:
            atom->covrad = 160/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_LU;
            atom->number = ATOM_LU;
            strcpy(atom->symbol, "Lu");
            atom->coreel = 54;
            atom->mass = 174.967;
            break;
        case ATOM_HF:
            atom->covrad = 150/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_HF;
            atom->number = ATOM_HF;
            strcpy(atom->symbol, "Hf");
            atom->coreel = 54;
            atom->mass = 178.49;
            break;
        case ATOM_TA:
            atom->covrad = 138/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_TA;
            atom->number = ATOM_TA;
            strcpy(atom->symbol, "Ta");
            atom->coreel = 54;
            atom->mass = 180.9479;
            break;
        case ATOM_W:
            atom->covrad = 146/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_W;
            atom->number = ATOM_W;
            strcpy(atom->symbol, "W");
            atom->coreel = 54;
            atom->mass = 183.84;
            break;
        case ATOM_RE:
            atom->covrad = 159/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_RE;
            atom->number = ATOM_RE;
            strcpy(atom->symbol, "Re");
            atom->coreel = 54;
            atom->mass = 186.207;
            break;
        case ATOM_OS:
            atom->covrad = 128/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_OS;
            atom->number = ATOM_OS;
            strcpy(atom->symbol, "Os");
            atom->coreel = 54;
            atom->mass = 190.23;
            break;
        case ATOM_IR:
            atom->covrad = 137/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_IR;
            atom->number = ATOM_IR;
            strcpy(atom->symbol, "Ir");
            atom->coreel = 54;
            atom->mass = 192.217;
            break;
        case ATOM_PT:
            atom->covrad = 128/100./BOHR;
            atom->rvdw = 175./100./BOHR;
            atom->charge = ATOM_PT;
            atom->number = ATOM_PT;
            strcpy(atom->symbol, "Pt");
            atom->coreel = 54;
            atom->mass = 195.078;
            break;
        case ATOM_AU:
            atom->covrad = 144/100./BOHR;
            atom->rvdw = 166./100./BOHR;
            atom->charge = ATOM_AU;
            atom->number = ATOM_AU;
            strcpy(atom->symbol, "Au");
            atom->coreel = 54;
            atom->mass = 196.9665;
            break;
        case ATOM_HG:
            atom->covrad = 149/100./BOHR;
            atom->rvdw = 155./100./BOHR;
            atom->charge = ATOM_HG;
            atom->number = ATOM_HG;
            strcpy(atom->symbol, "Hg");
            atom->coreel = 54;
            atom->mass = 200.59;
            break;
        case ATOM_TL:
            atom->covrad = 148/100./BOHR;
            atom->rvdw = 196./100./BOHR;
            atom->charge = ATOM_TL;
            atom->number = ATOM_TL;
            strcpy(atom->symbol, "Tl");
            atom->coreel = 54;
            atom->mass = 204.3833;
            break;
        case ATOM_PB:
            atom->covrad = 147/100./BOHR;
            atom->rvdw = 202./100./BOHR;
            atom->charge = ATOM_PB;
            atom->number = ATOM_PB;
            strcpy(atom->symbol, "Pb");
            atom->coreel = 54;
            atom->mass = 207.2;
            break;
        case ATOM_BI:
            atom->covrad = 146/100./BOHR;
            atom->rvdw = 207./100./BOHR;
            atom->charge = ATOM_BI;
            atom->number = ATOM_BI;
            strcpy(atom->symbol, "Bi");
            atom->coreel = 54;
            atom->mass = 208.9804;
            break;
        case ATOM_PO:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 197./100./BOHR;
            atom->charge = ATOM_PO;
            atom->number = ATOM_PO;
            strcpy(atom->symbol, "Po");
            atom->coreel = 54;
            atom->mass = 209;
            break;
        case ATOM_AT:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 202./100./BOHR;
            atom->charge = ATOM_AT;
            atom->number = ATOM_AT;
            strcpy(atom->symbol, "At");
            atom->coreel = 54;
            atom->mass = 210;
            break;
        case ATOM_RN:
            atom->covrad = 145/100./BOHR;
            atom->rvdw = 220./100./BOHR;
            atom->charge = ATOM_RN;
            atom->number = ATOM_RN;
            strcpy(atom->symbol, "Rn");
            atom->coreel = 54;
            atom->mass = 222;
            break;
        case ATOM_FR:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 348./100./BOHR;
            atom->charge = ATOM_FR;
            atom->number = ATOM_FR;
            strcpy(atom->symbol, "Fr");
            atom->coreel = 88;
            atom->mass = 223;
            break;
        case ATOM_RA:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 283./100./BOHR;
            atom->charge = ATOM_RA;
            atom->number = ATOM_RA;
            strcpy(atom->symbol, "Ra");
            atom->coreel = 88;
            atom->mass = 226;
            break;
        case ATOM_AC:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_AC;
            atom->number = ATOM_AC;
            strcpy(atom->symbol, "Ac");
            atom->coreel = 88;
            atom->mass = 227;
            break;
        case ATOM_TH:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_TH;
            atom->number = ATOM_TH;
            strcpy(atom->symbol, "Th");
            atom->coreel = 88;
            atom->mass = 232.0381;
            break;
        case ATOM_PA:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_PA;
            atom->number = ATOM_PA;
            strcpy(atom->symbol, "Pa");
            atom->coreel = 88;
            atom->mass = 231.0359;
            break;
        case ATOM_U:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 186./100./BOHR;
            atom->charge = ATOM_U;
            atom->number = ATOM_U;
            strcpy(atom->symbol, "U");
            atom->coreel = 88;
            atom->mass = 238.0289;
            break;
        case ATOM_NP:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_NP;
            atom->number = ATOM_NP;
            strcpy(atom->symbol, "Np");
            atom->coreel = 88;
            atom->mass = 237;
            break;
        case ATOM_PU:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_PU;
            atom->number = ATOM_PU;
            strcpy(atom->symbol, "Pu");
            atom->coreel = 88;
            atom->mass = 244;
            break;
        case ATOM_AM:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_AM;
            atom->number = ATOM_AM;
            strcpy(atom->symbol, "Am");
            atom->coreel = 88;
            atom->mass = 243;
            break;
        case ATOM_CM:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_CM;
            atom->number = ATOM_CM;
            strcpy(atom->symbol, "Cm");
            atom->coreel = 88;
            atom->mass = 247;
            break;
        case ATOM_BK:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_BK;
            atom->number = ATOM_BK;
            strcpy(atom->symbol, "Bk");
            atom->coreel = 88;
            atom->mass = 247;
            break;
        case ATOM_CF:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_CF;
            atom->number = ATOM_CF;
            strcpy(atom->symbol, "Cf");
            atom->coreel = 88;
            atom->mass = 251;
            break;
        case ATOM_ES:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_ES;
            atom->number = ATOM_ES;
            strcpy(atom->symbol, "Es");
            atom->coreel = 88;
            atom->mass = 252;
            break;
        case ATOM_FM:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_FM;
            atom->number = ATOM_FM;
            strcpy(atom->symbol, "Fm");
            atom->coreel = 88;
            atom->mass = 257;
            break;
        case ATOM_MD:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_MD;
            atom->number = ATOM_MD;
            strcpy(atom->symbol, "Md");
            atom->coreel = 88;
            atom->mass = 258;
            break;
        case ATOM_NO:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_NO;
            atom->number = ATOM_NO;
            strcpy(atom->symbol, "No");
            atom->coreel = 88;
            atom->mass = 259;
            break;
        case ATOM_LR:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_LR;
            atom->number = ATOM_LR;
            strcpy(atom->symbol, "Lr");
            atom->coreel = 88;
            atom->mass = 262;
            break;
        case ATOM_RF:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_RF;
            atom->number = ATOM_RF;
            strcpy(atom->symbol, "Rf");
            atom->coreel = 88;
            atom->mass = 261;
            break;
        case ATOM_DB:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_DB;
            atom->number = ATOM_DB;
            strcpy(atom->symbol, "Db");
            atom->coreel = 88;
            atom->mass = 262;
            break;
        case ATOM_SG:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_SG;
            atom->number = ATOM_SG;
            strcpy(atom->symbol, "Sg");
            atom->coreel = 88;
            atom->mass = 266;
            break;
        case ATOM_BH:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_BH;
            atom->number = ATOM_BH;
            strcpy(atom->symbol, "Bh");
            atom->coreel = 88;
            atom->mass = 264;
            break;
        case ATOM_HS:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_HS;
            atom->number = ATOM_HS;
            strcpy(atom->symbol, "Hs");
            atom->coreel = 88;
            atom->mass = 277;
            break;
        case ATOM_MT:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_MT;
            atom->number = ATOM_MT;
            strcpy(atom->symbol, "Mt");
            atom->coreel = 88;
            atom->mass = 268;
            break;
        case ATOM_DS:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_DS;
            atom->number = ATOM_DS;
            strcpy(atom->symbol, "Ds");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_RG:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_RG;
            atom->number = ATOM_RG;
            strcpy(atom->symbol, "Rg");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_CN:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_CN;
            atom->number = ATOM_CN;
            strcpy(atom->symbol, "Cn");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_UUT:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_UUT;
            atom->number = ATOM_UUT;
            strcpy(atom->symbol, "Uut");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_FL:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_FL;
            atom->number = ATOM_FL;
            strcpy(atom->symbol, "Fl");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_UUP:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_UUP;
            atom->number = ATOM_UUP;
            strcpy(atom->symbol, "Uup");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_LV:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_LV;
            atom->number = ATOM_LV;
            strcpy(atom->symbol, "Lv");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_UUS:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_UUS;
            atom->number = ATOM_UUS;
            strcpy(atom->symbol, "Uus");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_UUO:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_UUO;
            atom->number = ATOM_UUO;
            strcpy(atom->symbol, "Uuo");
            atom->coreel = 88;
            atom->mass = UNAVAIL;
            break;
        case ATOM_DUMM:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ATOM_UUO;
            atom->number = ATOM_UUO;
            strcpy(atom->symbol, "Du");
            atom->coreel = 0;
            atom->mass = UNAVAIL;
            break;
        case ATOM_WANN:
            atom->covrad = 50/100./BOHR;
            atom->rvdw = 100./100./BOHR;
            atom->charge = ATOM_WANN;
            atom->number = ATOM_WANN;
            strcpy(atom->symbol, "X");
            atom->coreel = 0;
            atom->mass = ZERO;
            break;
        default:
            atom->covrad = UNAVAIL/100./BOHR;
            atom->rvdw = 0./100./BOHR;
            atom->charge = ZERO;
            atom->number = 0;
            strcpy(atom->symbol, "XXX");
            atom->coreel = 0;
            atom->mass = UNAVAIL;
            break;
            }
}

void guess_acceptable_bond_lengths( real ** bondlengths )
{
    int i, j;
    atom_t ahi, ahj;
    char symdex[6] = "index";
    char dummy[MAXSTRLEN];

    for ( i=1; i<LAST_ATOM; i++ )
    {
        for ( j=1; j<LAST_ATOM; j++ )
        {
            sprintf( &dummy[0], "%i", i );
            assign_atom_parameters( &symdex[0], &dummy[0], &ahi );

            sprintf( &dummy[0], "%i", j );
            assign_atom_parameters( &symdex[0], &dummy[0], &ahj );

            // printf("%i %i\n", ahi.number , ahj.number );
            // printf("%f %f\n", ahi.covrad, ahj.covrad );

            bondlengths[i][j] = ahi.covrad + ahj.covrad;
        }
    }

}

void copy_atom(atom_t * anew, atom_t * aold)
{
    anew->covrad = aold->covrad;
    anew->rvdw = aold->rvdw;
    anew->charge = aold->charge;
    anew->number = aold->number;
    anew->mass = aold->mass;
    strcpy(anew->symbol, aold->symbol);
    anew->coreel = aold->coreel;

    int i;

    for ( i=0; i<DIM; i++ )
        anew->coords[i] = aold->coords[i];
}
