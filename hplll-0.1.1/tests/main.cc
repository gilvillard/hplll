/* Integer matrix nullspace test file  

Created Dim  7 avr 2013 16:54:03 CEST
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
#include "lehmer.cc"
#include "matgen.h"
#include "relations.h"

#include "tools.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  

  matrix<FP_NR<mpfr_t> > A;   // Input matrix 
  ZZ_mat<mpz_t> C;
 
  int r=6; 
  int s=6; 
  int n=r*s+1;

  int setprec=7000;
  mpfr_set_default_prec(setprec);

  gen3r2s(A,n,r,s);

  int found;
  
  // Alpha must be less than prec by a factor of ||F|| for having alpha bits
  found=relation_f<long, double>(C, A, setprec, 60, 1000, 20);

  if (found ==1)  print2maple(C,n,1);
	      
  //relation_f<long, double>(C, F, 244, 60, 800, 40, FPLLL,0.99);
  
  // relation_f(ZZ_mat<mpz_t>& C, const matrix<FP_NR<mpfr_t> > A, long alpha,
  //		   long confidence_gap = 60, long shift = 200, long increment = 20, int lllmethod = FPLLL, double delta = 0.99);
  
  //found = relation_lift(C, A, setprec);
     
  //relation_lll<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > (C, A, setprec, 10, FPLLL);  

      
  return 0;
}
