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



/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  typedef FP_NR<mpfr_t>   RT;
  typedef Z_NR<mpz_t>  ZT;
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT;  // fpLLL  

  // ---------------------------------------------------------------------
  { 
  
    cout << "************************************************************************** " << endl; 
    int n=10;
    int nbbits=20000;
    //int start=0,startsec=0;

    double delta=0.8;

    A.resize(n+1,n); 
    AT.resize(n,n+1);  
    AT.gen_intrel(nbbits);
    transpose(A,AT);


    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,TRANSFORM,DEF_REDUCTION);
    B.hlll(delta);

    ZZ_mat<mpz_t> U;
    U.resize(n,n);

    U=B.getU();

   

    ZZ_mat<mpz_t> V;
    V.resize(n,n);

    NTL_inv(V,U);

  } 

 
  return 0;
}
