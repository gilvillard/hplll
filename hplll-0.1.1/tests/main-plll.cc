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
#include "plll.h"

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
    int d=100;
    int n;
    int nbbits=400;
    
    int i,j;

    int start,startsec;


    //n=d+1;  A.resize(n,d);  AT.resize(d,n); AT.gen_intrel(nbbits);
    n=d; A.resize(n,d);  AT.resize(d,n);  AT.gen_uniform(nbbits);
    transpose(A,AT);
    for (i=1; i<n; i++)
      for (j=0; j<d; j++) 
	if (i==j) A(i,j)=1; else A(i,j)=0; 

    /*A.resize(d,d);
    AT.resize(d,d);
    AT.gen_ajtai(2.2);
    transpose(A,AT);
    print2maple(A,d,d);*/

    mpfr_set_default_prec(max(nbbits,2*d)+max(20,nbbits/10));
    //mpfr_set_default_prec(d+max(10,nbbits/10));

    PLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A);
    //PLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A);

    //print2maple(B.getbase(),n,d);

    start=utime();
    startsec=utimesec();

    B.hlll(0.99);
    
    start=utime()-start;
    startsec=utimesec()-startsec;

    /*int start;
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A,NO_TRANSFORM,DEF_REDUCTION);
    start=utime();
    C.hlll(0.99);
    start=utime()-start;
    cout << "   LLL " << start/1000 << " ms" << endl;*/
 
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T1.isreduced(0.9);

    cout << "   bits = " << nbbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   nblov plll " << B.nblov  << endl;
    cout << "   time plll: " << start/1000 << " ms" << endl;
    cout << "   time plll: " << startsec << " s" << endl;

    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A,NO_TRANSFORM,DEF_REDUCTION);

    start=utime();
    startsec=utimesec();

    C.hlll(0.99);

    start=utime()-start;
    startsec=utimesec()-startsec;


    cout << "   nblov hlll " << C.nblov  << endl;
    cout << "   time hlll: " << start/1000 << " ms" << endl;
    cout << "   time hlll: " << startsec << " s" << endl;

  } 

 
  return 0;
}
