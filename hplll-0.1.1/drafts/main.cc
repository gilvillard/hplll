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

using namespace hplll; 

/* ***********************************************

        Wrapper for using slll  

   ********************************************** */


  
template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
lll_wrap(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, double delta, int reduction_method=0) { 

  Timer time,ttot,tinit,tsize,th,thlll;
  time.clear();
  ttot.clear();
  tinit.clear();
  tsize.clear();
  th.clear();
  thlll.clear();
  
  int i,j,k;

  int n=A.getRows();
  int d=A.getCols();

  ZZ_mat<ZT> B;
    
  int d1=min(d,20);

  // Initial reduction
  // -----------------
  
  B.resize(n,d1);

  for (i=0; i<n; i++)
    for (j=0; j<d1; j++)
      B.Set(i,j,A(i,j));

  verboseDepth= 1;
  Lattice<ZT, FT,  MatrixZT, MatrixFT>  L(B,NO_TRANSFORM,reduction_method);

  time.start();

  L.hlll(delta);
  
  time.stop();
  tinit+=time;
  
  C=L.getbase();

  cout << endl << "Initial  time: " << tinit << endl << endl;

  ttot+=tinit;
  
  // Rest of the reductions, one dimension at a time
  // -----------------------------------------------
  
  for (k = d1+1; k<=d; k++) {    // Real dims (not indices) 

    time.start();
   
    B.resize(n,k);

    for (i=0; i<n; i++)
      for (j=0; j<k-1; j++)
	B.Set(i,j,C(i,j));

    for (i=0; i<n; i++)
      B.Set(i,k-1,A(i,k-1));

    if (verboseDepth > 0) 
      cout << "Discovering+ vector " << k  << "/" << d << endl;
     
	    
    

    // Size reduction of the last column
    // ---------------------------------
    
    Lattice<ZT, FT,  MatrixZT, MatrixFT>  LR(B,NO_TRANSFORM,reduction_method);

    LR.householder_r(0);
    LR.householder_v(0);
    
    for (i=1; i<k-1; i++) {
       LR.householder_r(i);
       LR.householder_v(i);
    }

    if (reduction_method < 1) 
      LR.hsizereduce(k-1);
    else 
      LR.seysenreduce(k-1);
   
    verboseDepth=1;
    
    time.stop();
    ttot+=time;
    if (verboseDepth >0) 
      cout << "     Size reduction: " << time << endl;


    // Reduction with the last column 
    // ------------------------------

    SLattice<ZT, FT,  MatrixZT, MatrixFT>  L(LR.getbase(),4,NO_TRANSFORM,reduction_method);
     
    time.start();

    verboseDepth=0; 
    L.hlll(delta,4,4,1000000);

    verboseDepth=1;
    time.stop();

    ttot+=time;
     
    if (verboseDepth >0) {
      cout << "     Phase+: " << time << endl;
      cout << "     Total: " << ttot << endl;
    }
   
    
    
    C=L.getbase();

    
    // Pour comparaison avec hplll
    // ---------------------------
    // Lattice<ZT, FT,  MatrixZT, MatrixFT>  LH(B,NO_TRANSFORM,reduction_method);

    // verboseDepth=0; 
    // th.start();
    
    // LH.hlll(delta);

    // th.stop();
    // thlll+=th;
    // cout << endl << "     hlll: " << th << endl << endl;
    // verboseDepth=1; 
    
  }

  cout << endl; 
  cout << "Initial reduction: " << tinit << endl;
  cout << "Segment reduction: " << ttot << endl;

  //cout << endl << endl << "lllw: " << tinit+ttot << endl; 

  //cout << endl << "hlll: " << tinit+thlll << endl << endl;
  
  return 0;


  } 

    
/* ***********************************************

          MAIN   

   ********************************************** */



int main(int argc, char *argv[])  {
  
  FP_NR<dd_real> tt;
  tt=1.0;
  
  ZZ_mat<mpz_t> A; // For hpLLL 

  ZZ_mat<mpz_t> AT;  // fpLLL
 
  // ---------------------------------------------------------------------
  { 
  
  int d=8;
  int n;
  
  //int nbbits=8;
  //double alpha;
  //int output;
  
  double delta = 0.99;

  //int m=1;
  
  //int lovmax=1000000;

  int S=4;  // Rajouter à commandline 

  int nbthreads=4;

  command_line_basis(A, n, d, delta, argc, argv);


  // With the wrapper
  // ----------------
  
  ZZ_mat<mpz_t> C; 
  

  Timer tw;
  tw.start();
     
  //lll_wrap<mpz_t, ldpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<long double, ldpe_t> > (C,A,delta,SEYSEN_REDUCTION);
  lll_wrap<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > (C,A,delta,DEF_REDUCTION);

  tw.stop();
  
  //cout << endl << "lllw: " << tw << endl; // cf bout de hlll aussi pour l'instant 

  
   // With hlll
   // ---------

  verboseDepth=0;
   
   // Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L(A,NO_TRANSFORM,DEF_REDUCTION);
   
   // Timer tl;
   // tl.start();
   // L.hlll(delta);
   // tl.stop();

  
   // cout << endl << "hlll: " << tl << endl;

   // Vertification
   // -------------

   cout << endl << endl;

   verboseDepth=0;
   
   Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Btest(C,NO_TRANSFORM,DEF_REDUCTION);
   Btest.isreduced(delta-0.1);

   // Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Ltest(L.getbase(),NO_TRANSFORM,DEF_REDUCTION);
   // Ltest.isreduced(delta-0.1);


   //DBG ratio 
   double t,u,v,w;

   ratio<mpz_t>(C,t,u,v,w);
   
   cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
   cout << ".. Average diagonal ratio: " << u << endl;
   cout << ".. Max diagonal ratio: " << v << endl;
   cout << ".. First vector quality: " << w << endl;

   // ZZ_mat<mpz_t> Ct;

   // Ct.resize(d,n);
   // transpose(Ct,C);
   
   // cout << Ct << endl;
   
  } 

  
    
  return 0;
}
