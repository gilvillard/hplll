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


#include "tools.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  

  // --------------------------------------------------------------------- 
  
  int transform=1;

  double rho = 1.0;


  int n,d,i;
  double delta;
  
  ZZ_mat<mpz_t> A;  // fpLLL  

  command_line_basis(A, n, d, delta, argc, argv); 
 
  ZZ_mat<mpz_t> AT;  // fpLLL  

  AT.resize(d,n);
  
  transpose(AT,A);
 

  //Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > P(A,NO_TRANSFORM,DEF_REDUCTION);
  int prelimin=utime();
  //lllReduction(AT, 0.4, 0.51, LM_FAST,FT_DEFAULT,0);
  prelimin=utime()-prelimin;

  //int prelimin=utime();
  //P.hlll(0.4);
  //prelimin=utime()-prelimin;

  //A = P.getbase();
  transpose(A,AT);
 
  //print2maple(A,n,d);

    // ---------------------------------------------------------
    // Nb bits to consider, mpfr lattice 

    int height;
    height = size_in_bits(A(2,2)) - size_in_bits(A(n-3,n-3)) +1;
    //height =  size_in_bits(A(0,0))-size_in_bits(A(1,0)); 
   
    int bits;
    bits =   (4* (n + height) +n);

    mpfr_set_default_prec(bits);

    
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A);


    for (i=0; i<d; i++) {
	
	B.hsizereduce(i);
	B.householder_v(i);
      }

    // !!! Cond détruit R, à refaire après 
    FP_NR<mpfr_t> cc;
    cc=B.lcond(ANY,height,CHECK);
    Z_NR<mpz_t> ccond;
    ccond.set_f(cc);

    long cond;
    cond = ccond.get_si();

    B.assign(A);

    for (i=0; i<d; i++) {
	
	B.hsizereduce(i);
	B.householder_v(i);
      }

    A = B.getbase();
    transpose(AT,A);


    // ---------------------------------------
    // Approximate lattice 
    // -------------------

   
    bits =   (long) (((double) cond +1) * rho);
    cout << " ** New bits = " << bits << endl; 

    matrix<Z_NR<mpz_t> > RZ;
    RZ.resize(d,d);
    
    set_f(RZ,B.getR(),bits);
   
    ZZ_mat<mpz_t> Rtrunc;
    Rtrunc.resize(d,d);

    set(Rtrunc,RZ);

    // ----------------------------------------
    //  Reductions 

    int start,startinter;


    // HLLL trunc 
    // ----------
    
    ZZ_mat<mpz_t> res2;
    res2.resize(n,d);
    
    int hllltime;
    int hlllprod; 
    int hlllss;

    if (transform ==1) {
      
      start = utime();
      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Btruncss(Rtrunc,NO_TRANSFORM,DEF_REDUCTION);
      Btruncss.hlll(delta);
      hlllss=utime()-start;

      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Btrunc(Rtrunc,TRANSFORM,DEF_REDUCTION);

      start=utime();
      
      Btrunc.hlll(delta);
    
      startinter=utime();
 
      matprod(res2,A,Btrunc.getU());
   
    
      hllltime=utime()-start;
      hlllprod=utime()-startinter;
      
    }
    else {

      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Btrunc(Rtrunc,NO_TRANSFORM,DEF_REDUCTION);
      
      start=utime();
      
      Btrunc.hlll(delta);
      
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
      T2.isreduced(delta-0.1);

    // FPLLL trunc  
    // -----------
    transpose(AT,A);
    
    ZZ_mat<mpz_t> V;
    V.resize(d,d);
    for (i=0; i<d; i++) 
      V(i,i)=1;

    ZZ_mat<mpz_t> RtruncT;
    RtruncT.resize(d,d);
    transpose(RtruncT,Rtrunc);

    int fplllss;
    start = utime();
    lllReduction(RtruncT, delta, 0.51, LM_FAST,FT_DEFAULT,0);
    fplllss=utime()-start;

    transpose(RtruncT,Rtrunc);

    int fpllltime;
    int fplllprod;

    ZZ_mat<mpz_t> res3;
    res3.resize(n,d);

    res2.resize(d,n);
 
    if (transform ==1) {
      start=utime();
 
      lllReduction(RtruncT, V, delta, 0.51, LM_FAST,FT_DEFAULT,0);
      
      startinter=utime();
      
      matprod(res2,V,AT);
      
      fpllltime=utime()-start;
      fplllprod=utime()-startinter;

      transpose(res3,res2);
    } 
    else {
      start=utime();
 
      res2=RtruncT;

      lllReduction(res2, delta, 0.51, LM_FAST,FT_DEFAULT,0);
      
      startinter=utime();
      
      fpllltime=utime()-start;
      fplllprod=utime()-startinter;

      transpose(res3,res2);
    }

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T3(res3,NO_TRANSFORM,DEF_REDUCTION);
    T3.isreduced(delta-0.1);


    // Direct reduction 
    // ----------------

    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > DA(A,NO_TRANSFORM,DEF_REDUCTION);


    start=utime();
    DA.hlll(delta);
    int dhllltime=utime()-start;
    
    start=utime();
    lllReduction(AT, delta, 0.51, LM_FAST,FT_DEFAULT,0);
    int dfpllltime=utime()-start;
    
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T4(DA.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T4.isreduced(delta-0.1);

    ZZ_mat<mpz_t> TT;
    TT.resize(n,d);
    transpose(TT,AT);

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T5(TT,NO_TRANSFORM,DEF_REDUCTION);
    T5.isreduced(delta-0.1);

    cout << " initial  total  size = " << maxbitsize(A) << endl; 
    cout << " truncated total size = " << maxbitsize(Rtrunc) << endl << endl;
    cout << " cond = " << cond <<  "    height = " << height << endl;
    cout << " bits = " << bits << endl; 
   
     
    cout << endl; 

    
    cout << "      hlllss: " << hlllss/1000 << " ms" << endl;
    cout << "   time hlll: " << hllltime/1000 << " ms" << endl;
    cout << "        prod: " << hlllprod/1000 << " ms" << endl;
    cout << "      fplllss: " << fplllss/1000 << " ms" << endl;
    cout << "   time fplll: " << fpllltime/1000 << " ms" << endl;
    cout << "         prod: " << fplllprod/1000 << " ms" << endl;
    cout << "   time direct hlll: " << dhllltime/1000 << " ms" << endl;
    cout << "   time direct fplll: " << dfpllltime/1000 << " ms" << endl;
    cout << "   preliminary time: " << prelimin/1000 << " ms" << endl;

    cout << endl; 
    cout << "Ratio fplll: " << ((double) dfpllltime)/((double) fpllltime) << endl;
    cout << "Ratio  hlll: " << ((double) dhllltime)/((double) hllltime) << endl;
cout << "Truncation ratio: " <<  ((double) maxbitsize(A))/((double) maxbitsize(Rtrunc)) << endl; 
    cout << "Time trunc ratio: " << ((double) dfpllltime)/((double) fplllss) << endl; 


  return 0;
}
