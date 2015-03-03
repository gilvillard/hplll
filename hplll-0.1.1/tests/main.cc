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
 
  int r=8; 
  int s=9; 
  int n=r*s+1;

  int setprec=3200;
  mpfr_set_default_prec(setprec);

  gen3r2s(A,n,r,s);

  print2maple(A,1,n);

  // ZZ_mat<mpz_t> L;
  // L.resize(1,n);

  // for (int j=0; j<n; j++) {

  //   A(0,j).mul_2si( A(0,j), setprec);
  //   L(0,j).set_f(A(0,j));
  // }

  int found; 
  int start;
  start=utime();
  
  //found=relation_lift<mpz_t, double, matrix<FP_NR<double> > > (C, L, 3800, 0.99);
  found = relation_d(C, A, setprec);
  
  start=utime()-start;

  cout << "   dimension = " << n  << endl;
  cout << "   time relation: " << start/1000 << " ms" << endl;
  
  if (found ==1) print2maple(C,1,n);

  
    
  return 0;
}
