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
 
  int r=9; 
  int s=9; 
  int n=r*s+1;

  int setprec=5000;
  mpfr_set_default_prec(setprec);

  gen3r2s(A,n,r,s);

  ZZ_mat<mpz_t> L;
  L.resize(1,n);

  FP_NR<mpfr_t> t;
  
  for (int j=0; j<n; j++) {
    t.mul_2si( A(0,j), setprec);
    L(0,j).set_f(t);
  }

   int start=utime();
   
   relation_lift<long, double>(C, A, setprec, 40, FPLLL);
    
   // relation_lift_d_z<long, double> (C, L,  setprec, 100, 0.99, HLLL);
    
   start=utime()-start;
   
   
   cout << "   time internal: " << start/1000 << " ms" << endl;
      
    return 0; 

  
  //found = relation_lift(C, A, setprec);

  //found=relation_lll<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > (C, A, setprec, 10, FPLLL);  

      
  return 0;
}
