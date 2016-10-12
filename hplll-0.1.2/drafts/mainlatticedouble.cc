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
  

  // ---------------------------------------------------------------------
  
  int d=8;
  int nbbits=8;
  
  double delta = 0.75;


  PARSE_MAIN_ARGS {
    MATCH_MAIN_ARGID("-d",d);
    MATCH_MAIN_ARGID("-bits",nbbits);
    SYNTAX();
  }



  int start; 
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT;  // fpLLL  

  A.resize(d+1,d);
  AT.resize(d,d+1);  

  AT.gen_intrel(nbbits);
  transpose(A,AT);


  //---------
  
  Lattice<mpz_t, dpe_t,  matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);
   
  start=utime();
  
  B.hlll(delta);
  
  start=utime()-start;
  
  cout << endl; 
  cout << "   dimension = " << d  << endl;
  cout << "   time hlll: " << start/1000 << " ms" << endl;
    
  Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Btest(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);
  Btest.isreduced(delta-0.1);

  //-----------
  
  ZZ_mat<double> Ap;
  Ap.resize(d+1,d);
  
  for (int i=0; i<d+1; i++)
    for (int j=0; j<d ; j++) 
      Ap(i,j)=(double) mpz_get_si(A(i,j).getData());


  Lattice<double, dpe_t,  matrix<Z_NR<double> >, MatrixPE<double, dpe_t> > Bp(Ap,NO_TRANSFORM,DEF_REDUCTION);
  
  
  start=utime();
  
  Bp.hlll(delta);
  
  start=utime()-start;
  
  cout << endl; 
  cout << "   dimension = " << d  << endl;
  cout << "   time long: " << start/1000 << " ms" << endl;


  print2maple(B.getbase(),d+1,d);
  print2maple(Bp.getbase(),d+1,d);
  
  // Ap=Bp.getbase();
  
  // ZZ_mat<mpz_t> C;
  // C.resize(d+1,d);

  // for (int i=0; i<d+1; i++)
  //   for (int j=0; j<d ; j++) 
  //     C(i,j)=Ap(i,j).getData();
  
  // Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Bptest(C,NO_TRANSFORM,DEF_REDUCTION);
  // Bptest.isreduced(delta-0.1);
  
  
  return 0;
  }
