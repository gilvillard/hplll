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
#include "matgen.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A0,A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------

  int n,d;
  double delta;

  command_line_basis(A0, n, d, delta, argc, argv); 

  A.resize(n,d);
  AT.resize(d,n);
  transpose(AT,A);

    int start,startsec;

    Timer time;

    int status;
    
    cout << "--------------  HLLL" << endl << endl; 
    start=utime();
    startsec=utimesec();
   
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A0,NO_TRANSFORM);

    time.start();
    //status=B.hlll(0.7);
    status=B.hlll(delta);
    time.stop();

    start=utime()-start;
    startsec=utimesec()-startsec;
  
    
    cout << "   dimension = " << d  << endl;
    cout << "   time A: " << start/1000 << " ms" << endl;
    time.print(cout);
    
    
    if (status ==0) {
      Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,NO_LONG);
      T1.isreduced(delta-0.1);
      }
    cout << endl; 

    cout << "--------------  FPLLL WRAPPER" << endl << endl; 
    transpose(AT,A0);

    start=utime();
    startsec=utimesec();
    time.start();
    lllReduction(AT, delta, 0.501, LM_WRAPPER);
    time.stop();
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    
    cout << "   dimension = " << d  << endl;
    cout << "   time B: " << start/1000 << " ms" << endl;
    cout << "   time B: " << time << endl;

    transpose(A,AT);
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(A,NO_TRANSFORM,DEF_REDUCTION);
    T2.isreduced(delta-0.1);

   

  return 0;
}
