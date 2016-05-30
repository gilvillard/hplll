/* 

Created Dim  7 avr 2013 16:54:03 CEST
Copyright (C) 2013-2016      Gilles Villard 

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


#include "hplll.h"

#include "wrappers.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A;

  ZZ_mat<mpz_t> C;
  
  ZZ_mat<mpz_t> AT;
 
  // ---------------------------------------------------------------------

  int n,d;
  double delta;

  command_line_basis(A, n, d, delta, argc, argv); 

 
  Timer th,tf;
  
  // HLLL ------------------------------------------
   
  
  //Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);

  // verboseDepth = 1;
  // th.start();
  // status=B.hlll(delta);
  // th.stop();
  
  verboseDepth = 2;
  
  th=hlll<mpz_t>(C, A, 0.99, false, true);
    
  //hlll<__int128_t>(C, A, 0.99, 0,0); 
  //hlll<long>(C, A, 0.99, false, true); 
 
  cout << endl; 

  cout << "--------------  FPLLL WRAPPER" << endl << endl;
  
  AT.resize(d,n);
  
  transpose(AT,A);
  
  tf.start();
  
  lllReduction(AT, delta, 0.501, LM_FAST, FT_DEFAULT,0,LLL_VERBOSE);
  
  tf.stop();
  
  // transpose(A,AT);
  
  // Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A,NO_TRANSFORM,DEF_REDUCTION);
  // verboseDepth=0;
  //T.isreduced(delta-0.1);

  //  cout << "-----------------------" << endl;

   cout << "HLLL: " << th << endl;
  
   cout << "FPLLL :" << tf << endl;
  

  return 0;
}
