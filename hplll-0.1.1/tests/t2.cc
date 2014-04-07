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


#include "matgen.cc"
#include "relations.cc" 

#include "nullspace.cc" 


/* ***********************************************

          MAIN   

   ********************************************** */


int main(int argc, char *argv[])  {

  
  int d,n,bitsize;

  int start;

  d=12; 
  n=2*d;
  bitsize=4000;

  ZZ_mat<mpz_t> A;
  A.resize(d,n);
  A.gen_uniform(bitsize);

  int shift=17000; 

  for (int i=0; i<12 ; i++) { 

    ZZ_mat<mpz_t> C;
    
    start=utime();
    //nullspace_lll<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > (C, A, 10, 0, HLLL); 
    nullspace_lll<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > (C, A, shift, 0, HLLL);
    start = utime()-start;

    shift+=1000;

    cout << "    Shift = " << shift << endl; 
    cout << "    Dimension: " << d << " x " << n << endl;
    cout <<  "    Nullspace bit size = " << maxbitsize(C) << endl; 
    cout <<  "    Time: " << start/1000 << " ms" << endl;
    cout << " ----------------------------------------------- " << endl; 
  } 

  return 0;
}
