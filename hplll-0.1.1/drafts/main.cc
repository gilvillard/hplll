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


#include "hplll.h"

#include "../src/lehmer.cc"

using namespace hplll; 


/* ***********************************************

          MAIN   

   ********************************************** */



int main(int argc, char *argv[])  {

  typedef mpz_t  ZT;
  //typedef long ZT;

  ZZ_mat<ZT> A; // For hpLLL 

  // ---------------------------------------------------------------------
   
  int d=8;
  int n;
  
  
  double delta = 0.99;

  command_line_basis(A, n, d, delta, argc, argv);

  // Lehmer
  // ------
  
  ZZ_mat<ZT> C; 
  
  Timer tleh;
  tleh.start();

  verboseDepth=1;
  //lehmer_lll<long,dpe_t, MatrixPE<double, dpe_t> > (C, A, delta, 20);
  //lehmer_lll<__int128_t, double, matrix<FP_NR<double> > > (C, A, delta, 5, true);
  //lehmer_lll<mpz_t, double, matrix<FP_NR<double> > > (C, A, delta, 10, true);
  //lehmer_lll<mpz_t,dpe_t, MatrixPE<double, dpe_t> > (C, A, delta, 10, true);
  lehmer_lll<long, double, matrix<FP_NR<double> > > (C, A, delta, 5,true);
    
  tleh.stop();
  
  cout << endl << "Lehmer: " << tleh << endl;  

  
   // With hlll
   // ---------

  verboseDepth=0;
   
  Lattice<ZT, dpe_t, matrix<Z_NR<ZT> >, MatrixPE<double, dpe_t> > L(A,NO_TRANSFORM,DEF_REDUCTION);
   
  Timer tl;
  tl.start();
  //L.hlll(delta);
  tl.stop();

   cout << endl << "hlll: " << tl << endl;

   // Verification
   // -------------

   cout << endl << endl;

   verboseDepth=0;
   
   Lattice<ZT, mpfr_t,  matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > Btest(C,NO_TRANSFORM,DEF_REDUCTION);
   Btest.isreduced(delta-0.1);

   // Lattice<ZT, mpfr_t,  matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > Ltest(L.getbase(),NO_TRANSFORM,DEF_REDUCTION);
   // Ltest.isreduced(delta-0.1);

   

   //DBG ratio 
   // double t,u,v,w;

   // ratio<ZT>(C,t,u,v,w);
   
   // cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
   // cout << ".. Average diagonal ratio: " << u << endl;
   // cout << ".. Max diagonal ratio: " << v << endl;
   // cout << ".. First vector quality: " << w << endl;

   cout << "-----------------------" << endl;

   cout << "Lehmer LLL: " << tleh << endl;
   //tw.print(cout);
   cout << "HLLL :" << tl << endl;
   //tl.print(cout);  

    
  return 0;
}
