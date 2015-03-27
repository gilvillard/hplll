/* HJLS and PSLQ relation algorithms 

Created Jeu  7 mar 2013 15:01:41 CET
        Ven 27 mar 2015 14:18:38 CET

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


#ifndef HPLLL_RELATIONS_H
#define HPLLL_RELATIONS_H

#include "hlll.h" 
#include "matmixed.h"

#include "lehmer.cc"

#include "decomp.h"
#include "decompz.h"

namespace hplll { 

  // Calls underlying long-double routines
  
  template<class ZT, class FT> 
    int relation_f(ZZ_mat<mpz_t>& C, const matrix<FP_NR<mpfr_t> > A, long alpha, long conshift, int lllmethod=HLLL); 
    
  template<class ZT, class FT> int relation_lift_d_z(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A,  int alpha, int shift,
						     double delta, int lllmethod=HLLL);

  template<class ZT, class FT> int  
    detect_lift_d(ZZ_mat<ZT>& U, ZZ_mat<mpz_t> L_in, ZZ_mat<FT> A_in_f, int& new_def, int def,
		  int target_def, FP_NR<mpfr_t>& new_quot, int shift, double delta, int lllmethod=HLLL);
    
  template<class ZT, class FT, class MatrixFT> int  
    relation_lift_z(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int alpha=0, double delta=0.99);

  template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
    relation_lll(ZZ_mat<ZT>& C,  const matrix<FP_NR<mpfr_t> > L, long setprec, long shift, int lllmethod=HLLL);
  
  template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
    relation_lll_z(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int alpha, int shift=0, double delta=0.99, int lllmethod=HLLL);


} // end namespace hplll


#include "relations.cc"


#endif 
