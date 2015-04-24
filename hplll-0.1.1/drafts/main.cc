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
#include "slll.h"




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


    // Make the dimension divisible by K
    // ---------------------------------

   
    int i,j;
    
    if (d%K !=0) {

     
      ZZ_mat<mpz_t> B; // For hpLLL 
      B.resize(n+K-d%K,d+K-d%K); 

      Z_NR<mpz_t> amax;
      amax=0;

      for (i=0; i<d; i++) 
	if (A(i,i).cmp(amax) > 0) amax=A(i,i);
	
      
      for  (i=0; i<n; i++) 
	for (j=0; j<d; j++) 
	  B(i,j)=A(i,j);

      for  (i=0; i<K-d%K; i++)
	B(n+i,d+i)=amax;

      A.resize(n+K-d%K,d+K-d%K);
      set(A,B);

      AT.resize(d+K-d%K,n+K-d%K);
      transpose(AT,A);
    
      n+=K-d%K;
      d+=K-d%K;
      
    }
    else { 

      AT.resize(d,n);
      transpose(AT,A);
    }

        
    Timer time;

#ifdef _OPENMP
    OMPTimer ptime;  
#else 
    Timer ptime;
#endif 

    // TODO with long double and dpe_t
    SLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A,TRANSFORM,DEF_REDUCTION);

    ptime.start();

    // Régler la valeur de condbits 
    B.hlll(delta,40, K,lovmax);

    ptime.stop();

    // ZZ_mat<double> Uf;
    // Uf.resize(d,d);
    // Uf=B.getU();

    // FP_NR<double> tf;
    
    // ZZ_mat<mpz_t> U;
    // U.resize(d,d);
    // for (j=0; j<d; j++)
    //   for (i=0; i<d; i++) {
    // 	tf = Uf(i,j).getData();  // Pour long double ou autre, vérifier et passer par set_z ? 
    // 	U(i,j).set_f(tf);
    //   }

    // matprod_in(A,U);
	
    // Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(A,NO_TRANSFORM,DEF_REDUCTION);
    // T1.isreduced(delta-0.1);

    

    // cout << "   dimension = " << d  << endl;
    // cout << "   nblov plll " << B.nblov  << endl;
    // cout << "   time plll: " << ptime << endl;
    // cout << "   input bit size: " << maxbitsize(A) << endl;
    

    // transpose(A,AT);

    // Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A,TRANSFORM,DEF_REDUCTION);


    // time.start();

    // C.hlll(delta);

    // time.stop();

    // cout << endl; 
    // cout << "   nblov hlll " << C.nblov  << endl;
    // cout << "   time hlll: " << time << endl;

    // transpose(AT,A);

    // time.start();

    // lllReduction(AT, delta, 0.51, LM_WRAPPER,FT_DEFAULT,0);

    // transpose(A,AT);

    // Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(A,NO_TRANSFORM,DEF_REDUCTION);
    // T2.isreduced(delta-0.1);

    // time.stop();
    // cout << "   time fplll: " << time << endl;
    

    

  return 0;
}
