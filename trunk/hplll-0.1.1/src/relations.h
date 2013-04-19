/* HJLS and PSLQ relation algorithms 

Created Jeu  7 mar 2013 15:01:41 CET
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

#include "decomp.cc"

#ifndef HPLLL_RELATIONS_H
#define HPLLL_RELATIONS_H

// ********   ATTENTION   ******************************
//
// Template RT but tested for mpfr .....
//
/* *************************************************************

   COMPUTING INTEGER RELATIONS FOR A FLOATING POINT MATRIX 

   Method 1: HJLS, starting from the matrix d x n 
              using the identiy for calling decomp with dec > 0 

   Method 2: 1) compute an associated floating FGAS 
             2) go through the FGAS decomposition      

  
   d x n 
   Hence nullity (returned) is <= nbrel < n-d

   ************************************************************ */



template<class RT,  class FT, class MatrixFT>  int relations_hjls(ZZ_mat<mpz_t>& C, const matrix<FP_NR<RT> >& A, long nbrel, long setprec); 








#endif 
