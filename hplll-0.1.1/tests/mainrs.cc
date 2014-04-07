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


#include "matgen.h"
#include "relations.h" 

#include "nullspace.h" 

#include "matmixed.h" 

using namespace hplll;

/* ***********************************************

          MAIN   

   ********************************************** */



int main(int argc, char *argv[])  {
  
  
  typedef FP_NR<mpfr_t>   RT;
  typedef Z_NR<mpz_t>  ZT;
  
  // Real matrix 
  // -----------
  
  matrix<RT> F;

  int r=5; 
  int s=5;
 
  int n=r*s+1;
  
  // 6 6  1030        6 7   1400    7 7    1660     7  8  2260   8  8  2825 
  //   8   9   3488    9 9   4351    9 10   5382    10 10   6644

  int setprec=1000;
  mpfr_set_default_prec(setprec);

  gen3r2s(F,n,r,s);

  ZZ_mat<mpz_t> Id;
  
  Id.resize(n,n);
  for (int i=0; i<n; i++) 
    Id(i,i)=1;


  ZZ_mat<mpz_t> C;
 
  int start=utime();

  relations_lll<mpz_t, dpe_t, MatrixPE<double, dpe_t> > (C, F, setprec, 0, 200);

  //relations_hjls<mpfr_t, dpe_t, MatrixPE<double, dpe_t> >(C,F,1,setprec);

  //restarting_hjls<mpfr_t, double, matrix<FP_NR<double > > >(C,F,1,setprec);

  print2maple(C,n,1);

  start = utime()-start;

  cout << endl << "   Time: " << start/1000 << " ms" << endl;
  


  return 0;
}
