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
#include "plll.h"



/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT;  // fpLLL  

  // ---------------------------------------------------------------------
  
    cout << "************************************************************************** " << endl; 

    int n,d;
    double delta;
    
    int K=4; 

    command_line_basis(A, n, d, delta, argc, argv); 

    
    unsigned int lovmax = 4294967295;

    PARSE_MAIN_ARGS {
      MATCH_MAIN_ARGID("-K",K);
      MATCH_MAIN_ARGID("-lovmax",lovmax);
      }


    AT.resize(d,n);
    transpose(AT,A);

    Timer time;

    time.start();

#ifdef _OPENMP
#pragma omp parallel for 
#endif 
      for (int k=0; k<4; k++) {

	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A);
	B.hlll(delta);

      }

      
#ifdef _OPENMP
#pragma omp barrier
#endif 

      time.stop();
      cout << " Time: " << time << endl; 

    /* -----------
    int cond; 
  
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(A);

    cond = L.lcond(TRIANGULAR_PROPER);

    //cout << " cond = " << B.lcond(TRIANGULAR_PROPER) << endl; 
    //cout << " cond = " << B.lcond(ANY, DEFAULT_PREC) << endl;
    //cout << " cond = " << B.lcond(ANY, 10, CHECK) << endl; 
   

    mpfr_set_default_prec(cond);

    Timer time;
       
    int start,startsec;
    
    PLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A);

    start=utime();
    startsec=utimesec();

    time.start();

    B.hlll(delta,K,lovmax);

    time.stop();
    
    start=utime()-start;
    startsec=utimesec()-startsec;

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T1.isreduced(delta-0.1);

    
    cout << "   dimension = " << d  << endl;
    cout << "   nblov plll " << B.nblov  << endl;
    cout << "   time plll: " << start/1000 << " ms" << endl;
    cout << "   time plll: " << startsec << " s" << endl;
    cout << "   time plll: " << time  << endl;
    
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A,NO_TRANSFORM,DEF_REDUCTION);

    start=utime();
    startsec=utimesec();

    C.hlll(delta);

    start=utime()-start;
    startsec=utimesec()-startsec;


    cout << "   nblov hlll " << C.nblov  << endl;
    cout << "   time hlll: " << start/1000 << " ms" << endl;
    cout << "   time hlll: " << startsec << " s" << endl;

    cout << "K " << K << endl; 
    */ 

  return 0;
}
