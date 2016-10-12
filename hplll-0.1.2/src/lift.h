/* LIFT LLL  

Created Jeu 16 avr 2015 14:11:32 CEST 
        

Copyright (C) 2015     Gilles Villard 

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


#ifndef HPLLL_LIFT_H
#define HPLLL_LIFT_H

#include "hlll.h" 


namespace hplll { 

  // Calls underlying long-double routines
  
  
  template<class ZT, class FT>
    int lift_z(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A,  int alpha,
			  long confidence_gap = 60, long shift = 200, long increment = 20,
			  int lllmethod = FPLLL, double delta = 0.99);

  template<class ZT, class FT> int  
    lift_f_z(ZZ_mat<ZT>& U, ZZ_mat<mpz_t> L_in, ZZ_mat<FT> A_in_f, int& new_def, int def,
		    int target_def, FP_NR<mpfr_t>& new_quot,
		    long confidence_gap = 60, long shift = 200, long increment = 20,
		    int lllmethod = FPLLL, double delta = 0.99);


  
} // end namespace hplll


#include "lift.cc"


#endif 
