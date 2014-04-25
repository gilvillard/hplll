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


#include "tools.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  typedef FP_NR<mpfr_t>   RT;
  typedef Z_NR<mpz_t>  ZT;
  

  // --------------------------------------------------------------------- 
  
  int transform=1;

  double llldelta=0.75;
   
  double rho = 1.0;

    int sN, sX, h, delta;

    sN = 1024;
    sX = 324;
    h = 15;
    delta = 3;

    
    PARSE_MAIN_ARGS {
      MATCH_MAIN_ARGID("-sN",sN);
      MATCH_MAIN_ARGID("-sX",sX);
      MATCH_MAIN_ARGID("-h",h);
      MATCH_MAIN_ARGID("-r",rho);
      MATCH_MAIN_ARGID("-delta",llldelta);
      SYNTAX();
    }

    Z_NR<mpz_t> N,X,tz,one;

    one = 1;

    N.randb(sN);
    X.randb(sX);

    int n;

    n =delta *h;

    ZZ_mat<mpz_t> A; // For hpLLL 
    ZZ_mat<mpz_t> AT;  // fpLLL  

    AT.resize(n,n);
    A.resize(n,n);

    int i,j,k;

    tz=one;
    for (i=h-1; i>=0; i--) {
 
     for (k=0; k<delta; k++) 
	AT(i*delta+k,i*delta+k)=tz;
     
     tz.mul(tz,N);
    }

    tz=one;
    for (i=0; i<n; i++) {
 
      AT(i,i).mul(AT(i,i),tz);
      tz.mul(tz,X);
    }

    int s;

    for (i=0; i<n; i++) 
      for (j=i+1; j<n; j++) {

	s=size_in_bits(AT(i,i));
	AT(j,i).randb(s);

    }

    transpose(A,AT);
    //print2maple(A,n,n);
    
    // ---------------------------------------------------------
    // Nb bits to consider, mpfr lattice 

    int height;
    height = size_in_bits(A(delta-1,delta-1)) - size_in_bits(A(n-delta,n-delta)) +1;
   
   
    int bits;
    bits =   (4* (n + height) +n);

    mpfr_set_default_prec(bits);
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A);

    for (i=0; i<n; i++) {
	
	B.hsizereduce(i);
	B.householder_v(i);
      }

    // !!! Cond détruit R, à refaire après 
    FP_NR<mpfr_t> cc;
    cc=B.cond();
    Z_NR<mpz_t> ccond;
    ccond.set_f(cc);

    long cond;
    cond = ccond.get_si();

    B.assign(A);

    for (i=0; i<n; i++) {
	
	B.hsizereduce(i);
	B.householder_v(i);
      }

    A = B.getbase();

    // ---------------------------------------
    // Approximate lattice 
    // -------------------

   
    bits =  (long) (((double) cond +1) * rho);
    cout << " ** New bits = " << bits << endl; 

    matrix<Z_NR<mpz_t> > RZ;
    RZ.resize(n,n);
    
    set_f(RZ,B.getR(),bits);
   
    //print2maple(RZ,n,n);

    ZZ_mat<mpz_t> Rtrunc;
    Rtrunc.resize(n,n);

    set(Rtrunc,RZ);

    // ----------------------------------------
    //  Reductions 

    int start,startinter;


    // HLLL trunc 
    // ----------
    
    ZZ_mat<mpz_t> res2;
    res2.resize(n,n);
    
    int hllltime;
    int hlllprod; 

    if (transform ==1) {
      
      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Btrunc(Rtrunc,TRANSFORM,DEF_REDUCTION);

      start=utime();
      
      Btrunc.hlll(llldelta);
    
      startinter=utime();
 
      matprod(res2,A,Btrunc.getU());
   
      hllltime=utime()-start;
      hlllprod=utime()-startinter;
      
    }
    else {

      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Btrunc(Rtrunc,NO_TRANSFORM,DEF_REDUCTION);
      
      start=utime();
      
      Btrunc.hlll(llldelta);
      
      startinter=utime();
      
      res2=Btrunc.getbase();
      
      hllltime=utime()-start;
      hlllprod=utime()-startinter;
      
      //ZZ_mat<mpz_t> TT;
      //TT.resize(n,n);
      //NTL_inv(TT,A);
      //print2maple(TT,n,n);

    }
      
      Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(res2,NO_TRANSFORM,DEF_REDUCTION);
      T2.isreduced(llldelta-0.1);

    // FPLLL trunc  
    // -----------
    transpose(AT,A);
    
    ZZ_mat<mpz_t> V;
    V.resize(n,n);
    for (i=0; i<n; i++) 
      V(i,i)=1;

    ZZ_mat<mpz_t> RtruncT;
    RtruncT.resize(n,n);
    transpose(RtruncT,Rtrunc);

    int fpllltime;
    int fplllprod;

    ZZ_mat<mpz_t> res3;
    res3.resize(n,n);

    if (transform ==1) {
      start=utime();
 
      lllReduction(RtruncT, V, llldelta, 0.51, LM_WRAPPER,FT_DEFAULT,0);
      
      startinter=utime();
      
      matprod(res2,V,AT);
      
      fpllltime=utime()-start;
      fplllprod=utime()-startinter;

      transpose(res3,res2);
    } 
    else {
      start=utime();
 
      res2=RtruncT;

      lllReduction(res2, llldelta, 0.51, LM_WRAPPER,FT_DEFAULT,0);
      
      startinter=utime();
      
      fpllltime=utime()-start;
      fplllprod=utime()-startinter;

      transpose(res3,res2);
    }

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T3(res3,NO_TRANSFORM,DEF_REDUCTION);
    T3.isreduced(llldelta-0.1);


    // Direct reduction 
    // ----------------

    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > DA(A,NO_TRANSFORM,DEF_REDUCTION);


    start=utime();
    DA.hlll(llldelta);
    int dhllltime=utime()-start;
    
    start=utime();
    lllReduction(AT, llldelta, 0.51, LM_WRAPPER,FT_DEFAULT,0);
    int dfpllltime=utime()-start;


    cout << " initial  total  size = " << maxbitsize(A) << endl; 
    cout << " truncated total size = " << maxbitsize(Rtrunc) << endl << endl;
    cout << " cond = " << cond << endl;
    cout << " bits = " << bits << endl; 
    cout << " n = " << n << "    h = " << h << "    N = " << N << "    X = " << X  << "    height = " << height << endl; 
     
    cout << endl; 

    cout << "   time hlll: " << hllltime/1000 << " ms" << endl;
    cout << "        prod: " << hlllprod/1000 << " ms" << endl;
    cout << "   time fplll: " << fpllltime/1000 << " ms" << endl;
    cout << "         prod: " << fplllprod/1000 << " ms" << endl;
    cout << "   time direct hlll: " << dhllltime/1000 << " ms" << endl;
    cout << "   time direct fplll: " << dfpllltime/1000 << " ms" << endl;
  return 0;
}
