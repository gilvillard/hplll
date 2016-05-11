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
#include <iostream>

#include <gmp.h>
#include <mpfr.h>
#include <fplll.h>

#include "givtimer.C"

#include "nr_Z_l128.h"

#ifndef __CYGWIN__
#define HPLLL_WITH_LONG_DOUBLE
#endif

#ifdef FPLLL_WITH_LONG_DOUBLE 
#define HPLLL_WITH_LONG_DOUBLE  // Long double in hplll require long double in fplll
#include "ldpe.h"
#endif
#include "nr-ld.cpp"  // Should be cleaned and in part in the test above 


// By default, but does nothing at level 0  
#define HPLLL_VERBOSE

#ifdef HPLLL_VERBOSE
int verboseDepth = 0;
#define HPLLL_INFO(x,y) {if (verboseDepth >= 0) cout << x << y << endl;}
#else
#define HPLLL_INFO(x,y)
#endif 



namespace hplll { 

#define DEF_REDUCTION 0
#define SEYSEN_REDUCTION 1 

#define WITH_LONG 1
#define NO_LONG 0

#define TRANSFORM 1
#define INV_T_TRANSFORM 2 
#define NO_TRANSFORM 0
#define HLLL 1
#define FPLLL 0 
#define L1 2

//#define HPLLL_VERBOSE(x) {if cout << x << endl;}


using namespace std;

using namespace fplll;

using namespace Givaro;

} // end namespace hplll


#endif
