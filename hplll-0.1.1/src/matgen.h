/* Some matrix constructors  

Created Sam  6 avr 2013 17:42:48 CEST 
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

#include  "defs.h"
#include  "tools.h"
#include  "mat.h"
#include  "matpe.h"

#ifndef MATGEN_H
#define MATGEN_H

namespace hplll { 

  /* ***********************************************

     GENERATION   

     ********************************************** */

  /* Using the defaulft mpfr precision */
  /* A row vector  */
  
  // One real vector for a relation 
  template<class RT> int gen3r2s(matrix<FP_NR<RT> >& B, int n, int r, int s); 
  
  // Two real vector for a simultaneous relation 
  template<class RT> int gen3r2s7t5u(matrix<FP_NR<RT> >& B, int n, int r, int s, int t, int u);
  
  // When testing prec and Seysen
  int genalpha(ZZ_mat<mpz_t>& B, int n, double alpha);
  
  
  void command_line_basis(ZZ_mat<mpz_t>& A, int& n, int& d, double &delta, int argc, char *argv[]);
  
} // end namespace hplll


#include "matgen.cc"

#endif 
