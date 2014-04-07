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
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------
  { 
  
    cout << "************************************************************************** " << endl; 
    int d=160;
    int nbbits=200;
    int start,startsec;

    double delta=0.8;

    A.resize(d+1,d); 
    tmpmat.resize(d+1,d); 
    AT.resize(d,d+1);  
    AT.gen_intrel(nbbits);
    transpose(A,AT);

    cout << "--------------  SeysenLLL" << endl << endl; 
    start=utime();
    startsec=utimesec();
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,SEYSEN_REDUCTION);
    B.hlll(delta);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << nbbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time A: " << start/1000 << " ms" << endl;
    cout << "   time A: " << startsec << " s" << endl;

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T1.isreduced(delta);

    cout << endl; 

    cout << "--------------  HLLL" << endl << endl; 
   

    start=utime();
    startsec=utimesec();
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A,NO_TRANSFORM,DEF_REDUCTION);
    C.hlll(delta);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << nbbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time B: " << start/1000 << " ms" << endl;
    cout << "   time B: " << startsec << " s" << endl;


    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(C.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T2.isreduced(delta);


    /*
    cout << endl; 
  
    transpose(AT,A);
 
    start=utime();
    startsec=utimesec();
    lllReduction(AT, delta, 0.5, LM_FAST,FT_DEFAULT,0);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << nbbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time C: " << start/1000 << " ms" << endl;
    cout << "   time C: " << startsec << " s" << endl;
  
    transpose(tmpmat,AT);
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T3(tmpmat,NO_TRANSFORM,DEF_REDUCTION);
    T3.isreduced(delta);
    */

  } 

 
  return 0;
}
