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

#include "tools.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> C; // For Lehmer
  ZZ_mat<mpz_t> AT;  // fpLLL  

  // ---------------------------------------------------------------------
  { 
  
  int d=8;
  int nbbits=100;
  int shift = 0;
  int alpha = 0;
  double delta = 0.99;


    PARSE_MAIN_ARGS {
      MATCH_MAIN_ARGID("-d",d);
      MATCH_MAIN_ARGID("-bits",nbbits);
      MATCH_MAIN_ARGID("-shift",shift);
      MATCH_MAIN_ARGID("-alpha",alpha);
      MATCH_MAIN_ARGID("-delta",delta);
      SYNTAX();
    }



  int start; 

  A.resize(d+1,d); 
  AT.resize(d,d+1);  
  AT.gen_intrel(nbbits);
  transpose(A,AT);

  //print2maple(A,d+1,d);

  
  ZZ_mat<mpz_t> A_up;
  A_up.resize(1,d); 

  for (int j=0; j<d; j++)
    A_up(0,j)=A(0,j);
    
  start=utime();

  //lift_lll<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > > (C, A_up, shift, alpha, delta);
  lift_lll<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > (C, A_up, shift, alpha, delta);

  
  start=utime()-start;

  cout << endl; 
  cout << "   dimension = " << d  << endl;
  cout << "   time lehmer: " << start/1000 << " ms" << endl;

  
  Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Btest(C,NO_TRANSFORM,DEF_REDUCTION);
  Btest.isreduced(delta-0.1);

  // FPLLL
  // -----
  
  start=utime();
  lllReduction(AT, delta, 0.51, LM_FAST,FT_DEFAULT,0);
  start=utime()-start;
    
  cout << endl; 
  cout << "   bits = " << nbbits << endl;
  cout << "   dimension = " << d  << endl;
  cout << "   time fplll: " << start/1000 << " ms" << endl;
  
  //transpose(A,AT);
  //Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Ctest(A,NO_TRANSFORM,DEF_REDUCTION);
  //Ctest.isreduced(delta);
  

  } 

 
  return 0;
}
