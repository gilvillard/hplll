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
   
  int n=6;
  int nbbits = 4;
  double rho=1.0;
  int shift=0;

   PARSE_MAIN_ARGS {
     MATCH_MAIN_ARGID("-n",n);
     MATCH_MAIN_ARGID("-bits",nbbits);
     MATCH_MAIN_ARGID("-shift",shift);
     MATCH_MAIN_ARGID("-r",rho);
     MATCH_MAIN_ARGID("-delta",llldelta);
     
     SYNTAX();
    }

   int i,j;

  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT;  // fpLLL  
  
  AT.resize(n,n);
  A.resize(n,n);

 
  AT.gen_uniform(nbbits);
  for (i=0; i<n/2; i++)
    for (j=0; j<n; j++)
      AT(j,i).randb(nbbits+shift);
  transpose(A,AT);

  //print2maple(A,n,n);
    
    // ---------------------------------------------------------
    // Nb bits to consider, mpfr lattice 


    int bits;
    bits =   (3* (n + shift) +n);

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

    // height 
    // ------

    for (i=0; i<n; i++) {
	
      
      B.hsizereduce(i);
      B.householder_v(i);
    }

    matrix<FP_NR<mpfr_t> > RR;
    RR.resize(n,n);
    RR=B.getR();
    FP_NR<mpfr_t>  maxquo,tt;
    maxquo=0.0;

    for  (j=0; j<n; j++) 
      for (i=0; i<n; i++) {
	tt.div(RR(j,j),RR(i,i));
	tt.abs(tt);
	if (tt.cmp(maxquo) > 0) maxquo=tt;
      }

    maxquo.log(maxquo);
    FP_NR<mpfr_t> c2;
    c2=2.0;
    c2.log(c2);
    maxquo.div(maxquo,c2);
  
    Z_NR<long> height;
    height.set_f(maxquo);
    // End height 
    


    B.assign(A);

    for (i=0; i<n; i++) {
	
	B.hsizereduce(i);
	B.householder_v(i);
      }

    A = B.getbase();

    print2maple(A,n,n);

    // ---------------------------------------
    // Approximate lattice 
    // -------------------

    //bits =  (long) (((double) 2*cond) * rho + n);
   
    bits =  (long) (((double) cond +1) * rho);
    
    

    matrix<Z_NR<mpz_t> > RZ;
    RZ.resize(n,n);
    
    matrix<FP_NR<mpfr_t> > Afp;
    Afp.resize(n,n);
   
    set(RZ,A);

    for (i=0; i<n; i++) 
      Afp.setcol(i,RZ.getcol(i),0,n);

    set_f(RZ,Afp,bits);
   
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
    int hlllss;
    int sizeofU;

    if (transform ==1) {

      start = utime();
      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Btruncss(Rtrunc,NO_TRANSFORM,DEF_REDUCTION);
      //Lattice<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > > Btruncss(Rtrunc,NO_TRANSFORM,DEF_REDUCTION);
      Btruncss.hlll(llldelta);
      hlllss=utime()-start;

      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Btrunc(Rtrunc,TRANSFORM,DEF_REDUCTION);

      //Lattice<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > > Btrunc(Rtrunc,TRANSFORM,DEF_REDUCTION);

      start=utime();
      
      Btrunc.hlll(llldelta);
    
      startinter=utime();
 
      matprod(res2,A,Btrunc.getU());
      sizeofU=maxbitsize(Btrunc.getU());
   
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

    int fplllss;
    start = utime();
    lllReduction(RtruncT, llldelta, 0.51, LM_FAST,FT_DEFAULT,0);
    fplllss=utime()-start;

    transpose(RtruncT,Rtrunc);

    int fpllltime;
    int fplllprod;

    ZZ_mat<mpz_t> res3;
    res3.resize(n,n);

    if (transform ==1) {
      start=utime();
 
      lllReduction(RtruncT, V, llldelta, 0.51, LM_FAST,FT_DEFAULT,0);
      
      startinter=utime();
      
      matprod(res2,V,AT);
      
      fpllltime=utime()-start;
      fplllprod=utime()-startinter;

      transpose(res3,res2);
    } 
    else {
      start=utime();
 
      res2=RtruncT;

      lllReduction(res2, llldelta, 0.51, LM_FAST,FT_DEFAULT,0);
      
      startinter=utime();
      
      fpllltime=utime()-start;
      fplllprod=utime()-startinter;

      transpose(res3,res2);
    }

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T3(res3,NO_TRANSFORM,DEF_REDUCTION);
    T3.isreduced(llldelta-0.1);


    // Direct reduction 
    // ----------------

    transpose(AT,A);

    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > DA(A,NO_TRANSFORM,DEF_REDUCTION);


    start=utime();
    DA.hlll(llldelta);
    int dhllltime=utime()-start;
    
    start=utime();
    lllReduction(AT, llldelta, 0.51, LM_FAST,FT_DEFAULT,0);
    int dfpllltime=utime()-start;

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T4(DA.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T4.isreduced(llldelta-0.1);
    
    ZZ_mat<mpz_t> TT;
    TT.resize(n,n);
    transpose(TT,AT);
    
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T5(TT,NO_TRANSFORM,DEF_REDUCTION);
    T5.isreduced(llldelta-0.1);
    
    cout << " initial  total  size = " << maxbitsize(A) << endl; 
    cout << " truncated total size = " << maxbitsize(Rtrunc) << endl << endl;
    cout << " size of U = " << sizeofU << endl << endl;

    cout << " cond = " << cond <<   " height = " << height << endl;
    cout << " bits = " << bits << endl; 
    cout << " n = " << n  << endl; 
     
    cout << endl; 

    cout << "      hlllss: " << hlllss/1000 << " ms" << endl;
    cout << "   time hlll: " << hllltime/1000 << " ms" << endl;
    cout << "        prod: " << hlllprod/1000 << " ms" << endl;
    cout << "      fplllss: " << fplllss/1000 << " ms" << endl;
    cout << "   time fplll: " << fpllltime/1000 << " ms" << endl;
    cout << "         prod: " << fplllprod/1000 << " ms" << endl;
    cout << "   time direct hlll: " << dhllltime/1000 << " ms" << endl;
    cout << "   time direct fplll: " << dfpllltime/1000 << " ms" << endl;

    cout << endl; 
    cout << "Ratio fplll: " << ((double) dfpllltime)/((double) fpllltime) << endl;
    cout << "Ratio  hlll: " << ((double) dhllltime)/((double) hllltime) << endl;
    cout << "Truncation ratio: " <<  ((double) maxbitsize(A))/((double) maxbitsize(Rtrunc)) << endl; 
    cout << "Time trunc ratio: " << ((double) dfpllltime)/((double) fplllss) << endl; 
  return 0;
}
