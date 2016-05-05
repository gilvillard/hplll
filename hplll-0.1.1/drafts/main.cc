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
lll_wrap(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, double delta) { 

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
    
  int d1=min(d,10);

  // Initial reduction
  // -----------------
  
  B.resize(n,d1);

  for (i=0; i<n; i++)
    for (j=0; j<d1; j++)
      B.Set(i,j,A(i,j));

  verboseDepth= 1;
  Lattice<ZT, FT,  MatrixZT, MatrixFT>  L(B,NO_TRANSFORM,SEYSEN_REDUCTION);

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

    //Lattice<ZT, FT,  MatrixZT, MatrixFT>  L(B,NO_TRANSFORM,DEF_REDUCTION);

    if (verboseDepth > 0) {
      cout << "Discovering+ vector " << k  << "/" << d << endl;
      //cout << "     Phase: " << cpu_discovered << endl;
      //cout << "     Total: " << cpu_tot << endl;
    }
	    
    
    SLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t>  > L(B,4,NO_TRANSFORM,SEYSEN_REDUCTION);

    
   

    // Size reduction of the last column
    // ---------------------------------
    
    L.householder_r(0);
    L.householder_v(0);
    
    for (i=1; i<k-1; i++) {
       L.householder_r(i);
       L.householder_v(i);
    }

    
    L.seysenreduce(k-1);

    verboseDepth=1;
    
    time.stop();
    ttot+=time;
    if (verboseDepth >0) 
      cout << "     Size reduction: " << time << endl;
 
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
    // Lattice<ZT, FT,  MatrixZT, MatrixFT>  LH(B,NO_TRANSFORM,DEF_REDUCTION);

    // verboseDepth=0; 
    // th.start();
    
    // LH.hlll(delta);

    // th.stop();
    // thlll+=th;
    // cout << "    hlll: " << th << endl;
    
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
     
  lll_wrap<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > (C,A,delta);

  tw.stop();
  
  //cout << endl << "lllw: " << tw << endl; // cf bout de hlll aussi pour l'instant 

  
   // With hlll
   // ---------
   
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


   
  // // Découpage rectangle de A
  
 
  // ZZ_mat<mpz_t>  T;
  // T.resize(n,K);
  
  // for (i=0; i<n; i++)
  //   for (j=0; j<K; j++)
  //     T(i,j)=A(i,j);

  
  // // Réduction du premier rectangle
  
  // Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > LT(T,NO_TRANSFORM,DEF_REDUCTION);
 
  // Timer lt;
  // lt.start();
  // LT.hlll(delta);
  // lt.stop();
  // cout << endl << "d: " << d << "    " <<  "K: " << K << endl << endl; 
  // cout << "1er bout: " << lt << "   " <<  LT.nbswaps << endl;
  
  // // On complète de 1

  // T=LT.getbase();
  
  // ZZ_mat<mpz_t>  B;
  // B.resize(n,K+1);
  
  // for (i=0; i<n; i++)
  //   for (j=0; j<K; j++)
  //     B(i,j)=T(i,j);

  //  for (i=0; i<n; i++)
  //    B(i,j)=A(i,K);

    
  //  // Avec plll
  //  // ----------
   
  //  // PLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > LP(B,NO_TRANSFORM,DEF_REDUCTION);
   
  //  // Timer lp;
  //  // lp.start();
  //  // LP.hlll(delta);
  //  // lp.stop();
  //  // cout << endl << "plll: " << lp << "   " <<  LP.nbswaps << endl;
  //  // cout << endl << "compteur: " << LP.compteur <<  endl;

  //  //print2maple(B,n,K+1); 
	     
  //  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > TT(B,NO_TRANSFORM,DEF_REDUCTION);

  //  Timer lp0,lp;


  //  lp0.start();
  //  TT.householder();
  //  TT.hsizereduce(K);
  //  lp0.stop();
  
   
  //  S=4;
  //  nbthreads=4;


  //  lp.start();
   
  //  SLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t>  > LP(TT.getbase(),S,NO_TRANSFORM,DEF_REDUCTION);
 
  //  LP.hlll(0.99,S,nbthreads,lovmax);

  //  lp.stop();
  //  cout << endl << "hsize: " << lp0  << endl;
  //  cout << endl << "slll: " << lp  << endl;

  //  // Avec hlll
  //  // ----------
   
  //  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > LB(B,NO_TRANSFORM,DEF_REDUCTION);
   
  //  Timer lb;
  //  lb.start();
  //  LB.hlll(delta);
  //  lb.stop();

  
  //  cout << endl << "hplll: " << lb << "   " <<  LB.nbswaps << endl;
  //  cout << endl << "compteur: " << LB.compteur <<  endl;

  //  // Avec fplll
  //  // ----------

  //  ZZ_mat<mpz_t> BT;  
  //  BT.resize(K+1,n);
  //  transpose(BT,B);
   
  //  Timer fp;
  //  fp.start();
  //  lllReduction(BT, delta, 0.501, LM_WRAPPER);
  //  fp.stop();
  //  cout << endl << "fplll: " << fp << endl;

   
  //  cout << endl;
   
  //  Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Btest(LP.getbase(),NO_TRANSFORM,DEF_REDUCTION);
  //  Btest.isreduced(delta-0.1);

  } 

 
  return 0;
}
