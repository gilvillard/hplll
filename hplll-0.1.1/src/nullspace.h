/* Integer matrix nullspace  

Created Created Jeu  7 mar 2013 15:01:41 CET  
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

#include "hlll.h" 
#include "decomp.h"
#include "decompz.h"
#include "matmixed.h"
#include "l1.h"
#include "relations.h"

#ifndef HPLLL_NULLSPACE_H
#define HPLLL_NULLSPACE_H

namespace hplll { 

// Full rank is generally assumed 
// ******************************


// Direct integer FGAS decomposition 
// ---------------------------------

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
nullspace_direct_decomp(ZZ_mat<ZT>& C, ZZ_mat<ZT> A); 


// FGAS generation + decomposition (PSLQ spirit)  
// ---------------------------------------------

// Not templated w.r.t. other matrices for the moment at this upper level 
template<class RT, class FT, class MatrixFT> int  
fgasgen(matrix<FP_NR<RT> >& F,  ZZ_mat<mpz_t> A, const long prec); 

} // end namespace hplll


#include "nullspace.cc"



#endif 
