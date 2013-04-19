/* Global definitions

Created Mar  2 avr 2013 12:34:17 CEST
Copyright (C) 2013      Gilles Villard 

This file is part of the hplll Library 

The hplll Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The hplll Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the hplll Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */



#ifndef HPLLL_DEFS_H
#define HPLLL_DEFS_H

#include <cmath>
#include <sstream>

#include <gmp.h>
#include <mpfr.h>
#include <fplll.h>

#include "ldpe.h"
#include "nr-ld.cpp"


#define DEF_REDUCTION 0
#define SEYSEN_REDUCTION 1 
#define TRANSFORM 1
#define INV_T_TRANSFORM 2 
#define NO_TRANSFORM 0
#define HLLL 1
#define FPLLL 0 

//#define HPLLL_VERBOSE(x) {if cout << x << endl;}

/* #define FPLLL_BEGIN_NAMESPACE namespace fplll {
#define FPLLL_END_NAMESPACE }

FPLLL_BEGIN_NAMESPACE  */ 

using namespace std;

using namespace fplll;


//FPLLL_END_NAMESPACE

#endif
