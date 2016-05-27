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


/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  //typedef mpz_t  ZT;
  typedef long ZT;

  
  ZZ_mat<mpz_t> A0;
  ZZ_mat<ZT> A; // For hpLLL 
  ZZ_mat<ZT> AT;  // fpLLL  

  // ---------------------------------------------------------------------

  int n,d;
  double delta;

  command_line_basis(A0, n, d, delta, argc, argv); 

  // Attention en 128 bits, mpfr get_si pas autrement
  matrix_cast(A,A0);  // temporaire avec ci-dessous

  AT.resize(d,n);
  
  transpose(AT,A);

  

  Timer th,tf;

    int status;
    
    cout << "--------------  HLLL" << endl << endl; 
    

    //Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A0,NO_TRANSFORM,SEYSEN_REDUCTION);
    Lattice<long, long double, matrix<Z_NR<long> >, matrix<FP_NR<long double> > > B(A,NO_TRANSFORM);
 
    verboseDepth = 1;
    th.start();
    status=B.hlll(delta);
    th.stop();

    
    cout << "   dimension = " << d  << endl << endl;

    th.print(cout);

    //cout << B.getbase() << endl; 
      
    if (status ==0) {
      Lattice<ZT, mpfr_t, matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,DEF_REDUCTION,NO_LONG);
      verboseDepth = 0;
      T1.isreduced(delta-0.1);
      }
    cout << endl; 

    cout << "--------------  FPLLL WRAPPER" << endl << endl;
    
    transpose(AT,A);

   
    tf.start();

    //lllReduction(AT, delta, 0.501, LM_FAST, FT_DEFAULT,0,LLL_VERBOSE);

    tf.stop();
  
    
    cout << "   dimension = " << d  << endl << endl;
   
    tf.print(cout);
 
    transpose(A,AT);
    Lattice<ZT, mpfr_t, matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > T2(A,NO_TRANSFORM,DEF_REDUCTION);
    //T2.isreduced(delta-0.1);

   cout << "-----------------------" << endl;

   cout << "HLLL: " << th << endl;
   //tw.print(cout);
   cout << "FPLLL :" << tf << endl;
   //tl.print(cout);
   

  return 0;
}
