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

        Wrapper for using slll, with gap 

   ********************************************** */


// Reduced until real column gap_status 
  
template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
lll_wrap_gap(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int gap_position, double delta, int reduction_method=0) {

  int i,j;

  int gap_status;

  int n=A.getRows();
  int d=A.getCols();

  // First part of the basis
  // -----------------------
  
  ZZ_mat<ZT> B;

  B.resize(n,gap_position);

  for (i=0; i<n; i++)
    for (j=0; j<gap_position; j++)
      B.Set(i,j,A(i,j));
  
  // Reduction of the first part
  // ---------------------------
  
  SLattice<ZT, FT,  MatrixZT, MatrixFT>  LB(B,4,NO_TRANSFORM,reduction_method);

  gap_status=LB.hlll(delta,4,4,1000000);

  // Recursively
  // -----------
  if (gap_status >=2) {

    lll_wrap_gap<ZT, FT,  MatrixZT, MatrixFT>(B,LB.getbase(),gap_status,delta,reduction_method);

    
  }
  // For second part directly
  // ------------------------
  else {

    B=LB.getbase();
  }

  // Reduction of the second part
  // ----------------------------
  
  C.resize(n,d);

  for (j=0; j<gap_position; j++)
    for (i=0; i<n; i++)
      C.Set(i,j,B(i,j));

  for (j=gap_position; j<d; j++)
    for (i=0; i<n; i++)
      C.Set(i,j,A(i,j));
  
  Lattice<ZT, FT,  MatrixZT, MatrixFT>  LC(C,NO_TRANSFORM,reduction_method);
  
  LC.hlll(0.99);
  
  C=LC.getbase();

  return 0;
}





/* ***********************************************

        Wrapper for using slll  

   ********************************************** */


  
template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
lll_wrap(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int dthreshold, double delta, int reduction_method=0) { 

  verboseDepth-=1;
  
  OMPTimer time,ttot,tinit,tsize,th,thlll;
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

  int gap_status;
  
  int d1=min(d,dthreshold);

  // Initial reduction
  // -----------------
  
  B.resize(n,d1);

  for (i=0; i<n; i++)
    for (j=0; j<d1; j++)
      B.Set(i,j,A(i,j));

 
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
    
    if (verboseDepth >= 0) 
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
    
    time.stop();
    ttot+=time;

        
    if (verboseDepth >= 0) 
      cout << "     Size reduction: " << time << endl;

    
    time.start();

    
    // Reduction with the last column and recursive w.r.t. gap_status 
    // --------------------------------------------------------------

    SLattice<ZT, FT,  MatrixZT, MatrixFT>  L(LR.getbase(),4,NO_TRANSFORM,reduction_method);

    gap_status=L.hlll(delta,4,4,1000000);

    
    if (gap_status >=2) 
      lll_wrap_gap<ZT, FT,  MatrixZT, MatrixFT>(C,L.getbase(),gap_status,delta,reduction_method);

    else
      C=L.getbase();

    time.stop();

    ttot+=time;

    
    if (verboseDepth >=0) {
      cout << "     Phase+: " << time << endl;
      cout << "     Nblov: " << L.nblov << endl; 
      cout << "     Total: " << ttot << endl;
    }

   

    // Pour comparaison avec hplll
    // ---------------------------
    // Lattice<ZT, FT,  MatrixZT, MatrixFT>  LH(B,NO_TRANSFORM,reduction_method);

   
    // th.start();
    
    // LH.hlll(delta);

    // th.stop();
    // thlll+=th;
    // cout << endl << "     hlll: " << th << endl << endl;
   

  } // End loop on k: extra columns 
  

  cout << endl; 
  cout << "Initial reduction: " << tinit << endl;
  cout << "Segment reduction: " << ttot << endl;

  //cout << endl << endl << "lllw: " << tinit+ttot << endl; 

  //cout << endl << "hlll: " << tinit+thlll << endl << endl;

  verboseDepth+=1;
  return 0;

  } 

    
/* ***********************************************

          MAIN   

   ********************************************** */



int main(int argc, char *argv[])  {

  typedef mpz_t  ZT;
  //typedef long ZT;

  ZZ_mat<mpz_t> A0; // For hpLLL
  
  ZZ_mat<ZT> A; // For hpLLL 

  ZZ_mat<ZT> AT;  // fpLLL
 
  // ---------------------------------------------------------------------
   
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

  command_line_basis(A0, n, d, delta, argc, argv);

  // Attention en 128 bits, mpfr get_si pas autrement 
  matrix_cast(A,A0);  // temporaire avec ci-dessous

  // Random unimodular preconditioning for the 512 example
  // -----------------------------------------------------

  // ZZ_mat<ZT> Q;
  // Q.resize(d,d);

  // for (int i=0; i<d; i++)
  //   Q(i,i)=1;
  
  // for (int i=0; i<d; i++)
  //   for (int j=i+1; j<d ; j++) 
  //     Q(i,j).randb(4);

  // for (int i=0; i<d; i++)
  //   for (int j=0; j<i ; j++) 
  //     Q(i,j)=0;
 
  // matprod_in(A,Q);
  
  // for (int i=0; i<d; i++)
  //   Q(i,i)=1;
  
  // for (int i=0; i<d; i++)
  //   for (int j=i+1; j<d ; j++) 
  //     Q(i,j)=0;
  
  //  for (int i=0; i<d; i++)
  //   for (int j=0; j<i ; j++) 
  //     Q(i,j).randb(4);

  // matprod_in(A,Q);
  
  
  
  // With the wrapper
  // ----------------
  
  ZZ_mat<ZT> C; 


  Timer tw;
  tw.start();

  verboseDepth=2;

  
  lll_wrap<ZT, dpe_t, matrix<Z_NR<ZT> >, MatrixPE<double, dpe_t> > (C,A,10,delta,SEYSEN_REDUCTION);
  //lll_wrap<ZT, ldpe_t, matrix<Z_NR<ZT> >, MatrixPE<long double, ldpe_t> > (C,A,20,delta,SEYSEN_REDUCTION);
  //lll_wrap<ZT, dpe_t, matrix<Z_NR<ZT> >, MatrixPE<double, dpe_t> > (C,T,100,delta,SEYSEN_REDUCTION);

  tw.stop();
  
  //cout << endl << "lllw: " << tw << endl; // cf bout de hlll aussi pour l'instant 

  //cout << "Transposed result" << endl;

  //cout << transpose(C) << endl; 
  
   // With hlll
   // ---------

  verboseDepth=0;
   
   Lattice<ZT, dpe_t, matrix<Z_NR<ZT> >, MatrixPE<double, dpe_t> > L(A,NO_TRANSFORM,DEF_REDUCTION);
   
   Timer tl;
   tl.start();
   //L.hlll(delta);
   tl.stop();

  
   cout << endl << "hlll: " << tl << endl;

   // Verification
   // -------------

   cout << endl << endl;
  
   Lattice<ZT, mpfr_t,  matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > Btest(C,NO_TRANSFORM,DEF_REDUCTION);
   Btest.isreduced(delta-0.1);

   //Lattice<ZT, mpfr_t,  matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > Ltest(L.getbase(),NO_TRANSFORM,DEF_REDUCTION);
   //Ltest.isreduced(delta-0.1);


   //DBG ratio 
   double t,u,v,w;

   ratio<ZT>(C,t,u,v,w);
   
   cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
   cout << ".. Average diagonal ratio: " << u << endl;
   cout << ".. Max diagonal ratio: " << v << endl;
   cout << ".. First vector quality: " << w << endl;

   cout << "-----------------------" << endl;

   cout << "SLLL: " << tw << endl;
   //tw.print(cout);
   cout << "HLLL :" << tl << endl;
   //tl.print(cout);  

    
  return 0;
}
