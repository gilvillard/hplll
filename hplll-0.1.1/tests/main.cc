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

  int setprec=4000;
  mpfr_set_default_prec(setprec);

  gen3r2s(A,n,r,s);

  int found; 
  //int start;
  //start=utime();

  // !!!!! Quand directement en flottant
  
  //found=relation_lift<mpz_t, double, matrix<FP_NR<double> > > (C, L, 3800, 0.99);
  found = relation_lift(C, A, setprec);
  
  if (found == 1) print2maple(C,1,n);

  cout << endl << "----------" << endl;
 
  found=relation_lll<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > (C, A, setprec);  

  
  //if (found ==1) print2maple(C,1,n);

  // !!!!! Temporaire en passant pas une matrice entière

  // start=utime()-start;

  // cout << "   dimension = " << n  << endl;
  // cout << "   time relation: " << start/1000 << " ms" << endl;
  
  // if (found ==1) print2maple(C,n,1);

  // if (found > 1) print2maple(C,n,found);
    
  return 0;
}
