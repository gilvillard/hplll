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


    int cond; 
  
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(A);

    cond = L.lcond(TRIANGULAR_PROPER);

    //cout << " cond = " << B.lcond(TRIANGULAR_PROPER) << endl; 
    //cout << " cond = " << B.lcond(ANY, DEFAULT_PREC) << endl;
    //cout << " cond = " << B.lcond(ANY, 10, CHECK) << endl; 
   
    L.setprec(2*cond);

    // Truncation of the input lattice 
    // -------------------------------

    for (int j=0; j<d; j++) {
      
      L.hsizereduce(j);
      L.householder_v(j);
    }
    
    matrix<Z_NR<mpz_t> > RZ;
    RZ.resize(n,d);

    set_f(RZ,L.getR(),cond);
    set(A,RZ);
    
    mpfr_set_default_prec(cond);

    Timer time;

#ifdef _OPENMP
    OMPTimer ptime;  
#else 
    Timer ptime;
#endif 

    PLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A);

    ptime.start();

    B.hlll(delta,K,lovmax);

    ptime.stop();

         
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T1.isreduced(delta-0.1);

    
    cout << "   dimension = " << d  << endl;
    cout << "   nblov plll " << B.nblov  << endl;
    cout << "   time plll: " << ptime << endl;

    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A,NO_TRANSFORM,DEF_REDUCTION);


    time.start();

    C.hlll(delta);

    time.stop();

    cout << endl; 
    cout << "   nblov hlll " << C.nblov  << endl;
    cout << "   time hlll: " << time << endl;

    transpose(AT,A);

    time.start();

    lllReduction(AT, delta, 0.51, LM_WRAPPER,FT_DEFAULT,0);

    time.stop();
    cout << "   time fplll: " << time << endl;


    

  return 0;
}
