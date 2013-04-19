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

  d=2; 
  n=6;
  bitsize=8;


  ZZ_mat<mpz_t> A;
  A.resize(d,n);
  A.gen_uniform(bitsize);

  ZZ_mat<mpz_t> C;

  nullspace_lll<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > (C, A, 8, 0, FPLLL); 

  cout <<  "Nullspace bit size = " << maxbitsize(C) << endl; 
  
  return 0;
}
