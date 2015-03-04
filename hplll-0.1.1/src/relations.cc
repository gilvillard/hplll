/* HJLS and PSLQ relation algorithms 

Created Jeu  7 mar 2013 15:01:41 CET
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



#ifndef HPLLL_RELATIONS_CC
#define HPLLL_RELATIONS_CC

namespace hplll {

  /***********************************************************************************

      Restarting 
      Restricted to doubles for the moment 
      Calls LLL with elementary lifts 

  **************************************************************************************/ 
  
  int relation_lift(ZZ_mat<mpz_t>& C, const matrix<FP_NR<mpfr_t> >& A, long setprec) {

    mpfr_set_default_prec(setprec);

    int n=A.getCols();

    ZZ_mat<mpz_t> L;
    L.resize(1,n);

    FP_NR<mpfr_t> t;
  
    for (int j=0; j<n; j++) {
      t.mul_2si( A(0,j), setprec);
      L(0,j).set_f(t);
    }

    int found;

    int start=utime();
    found=relation_lift_z<mpz_t, double, matrix<FP_NR<double> > > (C, L, setprec, 0.99);
  
  
    start=utime()-start;

  
    cout << "   time internal: " << start/1000 << " ms" << endl;
  
    return found;
    
  } 

  /************************************************
     Companion integer subroutine for above one 
     Restricted to doubles for the moment 
  *************************************************/    

  template<class ZT, class FT, class MatrixFT> int  
  relation_lift_z(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int alpha=0, double delta=0.99) { 

    int m,d;
    int i,j;
  
    m=A.getRows();
    d=A.getCols();

    ZZ_mat<ZT> A_in;
    A_in.resize(m+d,d);

    for (i=0; i<m; i++)
      for (j=0; j<d; j++)
	A_in(i,j)=A(i,j);

    for (i=0; i<d; i++)
      A_in(m+i,i)=1;

    int bitsize = maxbitsize(A,0,m,d);
  
   
    // For assigning the exact upper lattice part at each step 
    ZZ_mat<mpz_t> L;
    L.resize(m,d);
   
    // For assigning the truncated basis at each step -- dpe_t ?? 
    ZZ_mat<double> T;
    T.resize(d+1,d);
    Z_NR<ZT> t;
   
    ZZ_mat<double> Uf;
    Uf.resize(d,d);

    ZZ_mat<long int> U;
    U.resize(d,d);

    for (i=0; i<d; i++)
      U(i,i)=1; 

    int def = -bitsize;

    int target_def = -bitsize + alpha;
   
    int new_def;

    int found;

    FP_NR<FT> rel_bound;
   
    // Main loop on the shifts
    // -----------------------
    while (def < target_def) {

      //int start;

      for (i=0; i<m; i++) 
	for (j=0; j<d; j++) 
	  L(i,j)=A_in(i,j);
    
      //current_shift += shift;

      //start=utime();
         
      for (i=0; i<d; i++)
	for (j=0; j<d ; j++) 
	  T(m+i,j).getData()=A_in(m+i,j).get_d(); // cf pb for assigning a double to Z_NR<double> 

      // Lattice<double, dpe_t,  matrix<Z_NR<double> >, MatrixPE<double, dpe_t> > Bp(T,TRANSFORM,DEF_REDUCTION,1);
      Lattice<double, FT,  matrix<Z_NR<double> >,  MatrixFT> Bp(T,TRANSFORM,DEF_REDUCTION,1);

      Bp.assignL(L);
     
      found = Bp.detect_lift(delta,def,target_def,new_def,maxbitsize(A_in,1,d,d),rel_bound);
      
      // start=utime()-start;
      // cout << "   time A: " << start/1000 << " ms" << endl << endl;
      // cout << " Def: " << new_def << endl;
      // start=utime();
     
     
      def=new_def;
     
      Uf=Bp.getU();

      for (i=0; i<d; i++)
	for (j=0; j<d ; j++) 
	  U(i,j).getData()=((long int) Uf(i,j).getData());
     
      matprod_in_si(A_in,U);

      if (found == 1) {
       
	C.resize(1,d);
	for (j=0; j<d; j++)
	  C(0,j)=A_in(m+j,0);
	return 1;
      }
     
      // start=utime()-start;
      // cout << "   time B: " << start/1000 << " ms" << endl << endl;
      // start=utime();

      // cout << "***** " << new_def << endl;
      // cout << "***** " << maxbitsize(A_in)+new_def << endl;
     
      // cout << endl << "  size of U: " << maxbitsize(U,0,d,d)  << endl;
     

    }

    // found = 0
    cout << "**** There might not be relations of norm less than " << rel_bound << endl; 
  
    return 0;
  

  } 



/* ***********************************************

   Lehmer like integer relations using LLL  

   d x n  input matrix 

   B m x n global lattice with the initial identity 

	  Returns the nullity as soon as >= n-d 
          Full rank required ? 

   alpha : bits, global shit for the upper part 
   lsigma : bits, elementary shift (step) size à la Lehmer 
   unique shift if lsigma=0;

   method : HLLL or FPLLL 

   ********************************************** */

/// !!!!!!!!!!!!!   FAIRE VERSION SANS TRONCATURE 
  
template<class ZT, class FT, class MatrixFT> int  
relations_lll(ZZ_mat<ZT>& C, matrix<FP_NR<mpfr_t> > F, int prec, int lalpha, int lsigma) { 

  
 mpfr_set_default_prec(prec); 

  int d,n,m;
  int i=0,j=0;

   d=F.getRows();
   n=F.getCols();
   m=d+n;

   int nullity=0;
   bool zeroedcol; 
   
   vector<int> nullcols;  // Zeroed columns in the upper part : 1, otherwise 0
   nullcols.resize(n); 

   // Copy in the temporary basis and initial shift 
   // ---------------------------------------------
   // ICI 
   //print2maple(F,d,n);

   MatrixRZ<matrix, FP_NR<mpfr_t>,Z_NR<ZT> > B;

   B.resize(d,n,n);


   for (i=0; i<d; i++) 
      for (j=0; j<n; j++) 
	B.setRT(i,j,F.get(i,j));

   Z_NR<ZT> one;
   one=1;
   for (i=0; i<n; i++) 
     B.setZT(i,i,one);

   B.shift(lalpha);
 
   // For the zero test 
   // -----------------
   
   FP_NR<mpfr_t> epsilon;
   epsilon=B.get_epsilon();
   cout << endl << "epsilon = " << epsilon << endl << endl; 

   epsilon=B.shift_epsilon(20);
   //epsilon=B.shift_epsilon(1000);  // Pour B4 

   cout << endl << "shifted epsilon = " << epsilon << endl << endl;

   FP_NR<mpfr_t> tmpabs; 

   // Unique shift 
   // ------------
   
   if (lsigma ==0) {  
     
     Lattice<ZT, FT, MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> >, MatrixFT> D(B,NO_TRANSFORM,DEF_REDUCTION);
       
     D.hlll(0.99);
       
       B=D.getmixedbase();

       B.shift(-lalpha);

       //print2maple(B,n+1,n);

       // Zero column test 
       // ----------------
       
       nullity=0;
   

       for (j=0; j<n; j++) {

	 nullcols[j]=0; 
	 zeroedcol=true;
	 for (i=0; i<d; i++) {
	   tmpabs.abs(B.getRT(i,j));
	   // ICI 
	   cout << tmpabs << "   vs   " << epsilon << endl; 
	   if (tmpabs.cmp(epsilon) >0 )  zeroedcol=false;  
	 }
	 if (zeroedcol) { 
	   nullity+=1;
	   nullcols[j]=1; 
	 } 
       }
     

     // Shift too small for being a unique shift 
     // ----------------------------------------
     
     if ((nullity==0) || (nullity==n)) { // attention dimensions 
       cout << endl << " *** Anomaly? Nullity is "  << nullity << ", too small shift or lack of precision?" << endl; 
       return -1; 
     }
     else {

       C.resize(n, nullity);
       
       int decj=0;
       for (j=0; j<n; j++) {  
	 if (nullcols[j]==1) {
	   for (i=0; i<n; i++) 
	     C(i,decj)=B.getZT(i,j);
	   decj++; 
	 }
       }

     } // Nullity ok 
     
     // ICI 
     cout << "  ****  Nb swaps: " <<  D.nblov << endl;  
     return nullity; 
     
   } // End unique shift 
   

   // Several shifts Lehmer like loop 
   // -------------------------------
   
   else {

     MatrixRZ<matrix, FP_NR<mpfr_t>,Z_NR<ZT> > BL;

     BL.resize(d,n,n);

     int max2up=0;
     max2up=B.maxbitsizeRT();

     Lattice<ZT, FT, MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> >, MatrixFT> Bt(B,TRANSFORM,DEF_REDUCTION);

     bool notfound=1;

     int current_shift;
     current_shift = -max2up;

     // Loop of elementary shifts 
     // -------------------------
    
     while (notfound) {
     //for (int K=0; K<1; K++) { 

       BL.set(B);

       current_shift+=lsigma;
       // ICI       
       cout << "current shift " << current_shift << endl; 
       
       BL.shift(current_shift);

       Bt.mixed_put(BL, lsigma+2*n); // Heuristic 
       
       Bt.hlll(0.9999);
       
       // Required multiplication update when BL has been truncated 
       matprod_in(B,Bt.getU());
       // if not truncated one should not multiply evrything        
      
       // Zero column test 
       // ----------------
       
       nullity=0;

       // ICI 
       //cout << endl << " ------------------------------- " << endl; 

       for (j=0; j<n; j++) {
	 nullcols[j]=0; 
	 zeroedcol=true;
	 for (i=0; i<d; i++) {
	   tmpabs.abs(B.getRT(i,j));
	   // ICI 
	   //cout << tmpabs << "   vs   " << epsilon << endl; 
	   if (tmpabs.cmp(epsilon) >0 )  zeroedcol=false;  
	 }
	 if (zeroedcol) { 
	   nullity+=1;
	   nullcols[j]=1; 
	 } 
       }

       if (nullity > 0) {

	 notfound=0;

	 C.resize(n, nullity);
	 
	 int decj=0;
	 for (j=0; j<n; j++) {  
	   if (nullcols[j]==1) {
	     for (i=0; i<n; i++) 
	       C(i,decj)=B.getZT(i,j);
	     decj++; 
	   }
	 }

	 // ICI 
	 //print2maple(C,n,nullity);
	 return nullity; 
       } 
       
     } // End while loop on the Lehmer shifts 

   } // End else Lehmer case 
   
   
   return nullity;
};


  
  
// ********   ATTENTION   ******************************
//
// Template RT but written for mpfr .....
//
/* *************************************************************

   COMPUTING INTEGER RELATIONS FOR A FLOATING POINT MATRIX 

   Method 1: HJLS, starting from the matrix d x n 
              using the identiy for calling decomp with dec > 0 

   Method 2: 1) compute an associated floating FGAS 
             2) go through the FGAS decomposition      

  
   d x n 
   Hence nullity (returned) is <= nbrel < n-d */

   

/***********************************************************************************

   HJLS 

   Starting from the matrix d x n using the identiy for calling decomp with dec > 0 

**************************************************************************************/ 


// Not templated w.r.t. other matrices for the moment at this upper level 
// Put the integer matrices 
template<class RT, class FT, class MatrixFT>  int relations_hjls(ZZ_mat<mpz_t>& C, const matrix<FP_NR<RT> >& A, long nbrel, long setprec) {

  int d,n;

  d=A.getRows();
  n=A.getCols();

  mpfr_set_default_prec(setprec); // Same as for decomp? 

  matrix<FP_NR<RT> > B;

  B.resize(n,n+d);

  int i,j;

  for (i=0; i<n; i++) {
    for (j=0; j<d; j++)
      B(i,j)=A(j,i);  // = is allowed, set is not needed? 
    B(i,d+i)=1;
  }



  Fgas<RT, mpz_t, FT,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, MatrixFT > F(B, NO_TRANSFORM, setprec,d);
     
  int flag_decomp;

 
  flag_decomp=F.decomp(1.1547005384, nbrel, LONG_MAX);
  //flag_decomp=F.decomp(1.4142135624, nbrel, LONG_MAX);

  if (flag_decomp==-1) {
    return -1;
  }


  // The nullspace C should be in the last n+d-rel columns of V=Transpose(U^(-1))
  // See analogous procedure in nullspace indirect (on other entry types) 
  // ******************************************************************************

  ZZ_mat<mpz_t> UT;
  UT.resize(n+d,n+d);
  transpose(UT,F.getV());

  C.resize(n,nbrel);

  for (i=0; i<n; i++)
    for (j=0; j<nbrel; j++) 
      C(i,j)=UT(n+d-nbrel+j,d+i); 

  // Check of the nullity 
  // Coul maybe not be used id decomp has succeeded?  
  // ***********************************************

  matrix<FP_NR<RT> > Cfp, ZM;

  Cfp.resize(n,nbrel);
  for (i=0; i<n; i++)
    for (j=0; j<nbrel; j++)
      set_z(Cfp(i,j),C(i,j));

  ZM.resize(d,nbrel);

  

  matprod(ZM,A,Cfp); // Should be epsilon-zero 
  

  /* Tuning of the epsilon = "zero" value for the test, about the floating point precision */
  /* ***************************************************************************************/

  FP_NR<RT> epsilon,tmp1;

  epsilon=F.epsilon;   // Should/could use a different epsilon? Yes if the precision of the input A is greater 


  /* Zero test    */
  /****************/

  int nullity=0;
  bool zeroedcol;

  
  for (j=0; j<nbrel; j++) {
       
    zeroedcol=true;

    for (i=0; i<d; i++) {
    
      tmp1.abs(ZM(i,j));  
      if (tmp1.cmp(epsilon) > 0)   
	zeroedcol=false;
    }
    if (zeroedcol) 
      nullity+=1;   
  }


  if (nullity != nbrel) {
    cout << " *******  ERROR: Input matrix non-full row rank or problem with the decomposition/precision" << endl; 
    cout << " *******          nullity = " << nullity << " <> " << nbrel << endl;  
    return -1;
  } 


  return(nbrel);

} // End body of FP relations HJLS  



   

/***********************************************************************************

   RESTARTING HJLS 

   Starting from the matrix d x n using the identiy for calling decomp with dec > 0 

**************************************************************************************/ 


// Not templated w.r.t. other matrices for the moment at this upper level 
// Put the integer matrices 
template<class RT, class FT, class MatrixFT>  int restarting_hjls(ZZ_mat<mpz_t>& C, const matrix<FP_NR<RT> >& A, long nbrel, long setprec) {

  int d,n;

  int nblov=0;

  d=A.getRows();
  n=A.getCols();

  int N;
  N=n+d;

  int tps=0; 
  int start; 


  FP_NR<RT> tr,acc; // Tmp  

  int i,j,k;
 
  mpfr_set_default_prec(setprec); // Same as for decomp?

  // Main real matrix 
  // ----------------
 
  // Main FGAS 
  matrix<FP_NR<RT> > B;
  B.resize(n,N);

  for (i=0; i<n; i++) {
    for (j=0; j<d; j++)
      B(i,j)=A(j,i);  // = is allowed, set is not needed? 
    B(i,d+i)=1;
  }

  matrix<FP_NR<RT> > Res;
  Res.resize(n,d);
 
  matrix<FP_NR<RT> > H;
  H.resize(n,N);

 
  // Float matrices 
  // --------------
  MatrixFT Bft;

  Bft.resize(n,N);

  MatrixFT Resft;

  Resft.resize(n,d);


  // Transform and inverse transpose 
  // -------------------------------

  // Global 
  ZZ_mat<mpz_t> U; 

  U.resize(N,N);
  for (i=0; i<N; i++)
    U(i,i)=1;

  // Local decomposition transformation 
  //ZZ_mat<mpz_t> FU; 
  ZZ_mat<long> FU;

  // Temporary FU in real format 
  matrix<FP_NR<RT> > FUr;
  FUr.resize(N,N);

  // Global inverse transpose transformation 
  ZZ_mat<mpz_t> V; 
  
  V.resize(N,N);
  for (i=0; i<N; i++)
    V(i,i)=1;
  
  // For later zero-epsilon test at high precision 
  FP_NR<RT> epsilon; 


  // Initial size reduction, especially for input cases
  // with very different magnitudes
  // And high precision FGAS, first loop preparation 

 

  Fgas<RT, mpz_t, RT,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, matrix<FP_NR<RT> > > G(B,NO_TRANSFORM,setprec);

  // ICI ce shift only for 1010 pb understanding 
  //G.shift_epsilon(220);
  epsilon=G.epsilon; 

  for (j=0; j<N; j++) 
    G.hsizereduce(j);

  matprod_in(V,G.getV());
  B=G.getfgas();
  H=G.getR();

  for (j=0; j<N; j++)
    Bft.setcol(j,H.getcol(j),0,n);
   


  // Residue update for next round 
  // -----------------------------

  for (i=0; i<n; i++)
    for (j=0; j<d; j++) {
      
      acc=0.0;
      for (k=0; k<n; k++) {
	set_z(tr,V(d+k,d+i));
	tr.mul(tr,A(j,k));  
	acc.add(acc,tr);
      }
      Res.set(i,j,acc);
    } 
  
  // Residue normalization  
  acc=0.0;
  for (i=0; i<n; i++)
    for (j=0; j<d; j++) {
      tr.abs(Res(i,j));
      if (tr.cmp(acc)>0) acc=tr;
    }
  
  for (j=0; j<d; j++)
    Res.div(j,0,Res.getcol(j,0),acc,n);
  
  for (j=0; j<d; j++)
    Resft.setcol(j,Res.getcol(j),0,n);
  

  // Main loop preparation 
  // ---------------------

  FP_NR<RT> r1;

  int flag_decomp;
  int nbloops=0;


  bool notfound=1;
  long internal_prec;
  internal_prec=(Resft.get(0,0)).getprec();


  // Low precision FGAS
  // ------------------

  Fgas<FT, long, FT,  MatrixFT, matrix<Z_NR<long> >, MatrixFT > F(Bft, Resft, TRANSFORM, internal_prec, d);

   
  F.shift_epsilon(16);
  F.shift_testU(4);
  // ICI que pour 1010 
  //F.shift_epsilon(12);
  //F.shift_testU(5);

  // ****************************************
  // Main loop on low precision decomposition
  // ****************************************


  while (notfound) {
  
    nbloops+=1;

    // FGAS decomposition  
    // ------------------

    F.set(Bft, Resft);
   
    start=utime();
 
    // ICI   restarting 
    flag_decomp=F.decomp(1.1547005384, nbrel, LONG_MAX);
    //flag_decomp=F.decomp(1.4142135624, nbrel, LONG_MAX);

    mpfr_set_default_prec(setprec); // Should be only an internal change above (but changes everything) 

    start = utime()-start;
    tps+=start;
    
    nblov+=F.nblov;

    
    
    // If relation found 
    // -----------------
    if (flag_decomp==0) {
    
      matprod_in_si(V,F.getV());
      notfound=0;

      // ICI 
      cout << endl << " WITH SMALL DECOMP ENDING " << endl << endl; 

      // ICI 
      print2maple(V,n+d,n+d);

    }
    else { // Relation not found, will continue 
           // ---------------------------------

      FU=F.getU();
      
      for (i=0; i<N; i++) 
	for (j=0; j<N; j++) {
	  set_z(tr,FU(i,j));
	  FUr.set(i,j,tr);
	}

     
      // Accumulate in B 
      // ---------------
      
      matprod_in(B,FUr);
      
      // V accumulation 
      // --------------
      
      matprod_in_si(V,F.getV());

      // and size-reduction for better numerical quality 
      // ----------------------------------------------- 

      G.set(B);

      for (j=0; j<N; j++) {
	G.hsizereduce(j);
      }

      matprod_in(V,G.getV());
      B=G.getfgas();

      H=G.getR();

      for (j=0; j<N; j++)
      Bft.setcol(j,H.getcol(j),0,n);

      
      // Residue update for next round 
      // -----------------------------

      for (i=0; i<n; i++)
	for (j=0; j<d; j++) {
	  
	  acc=0.0;
	  for (k=0; k<n; k++) {
	    set_z(tr,V(d+k,d+i));
	    tr.mul(tr,A(j,k));  
	    acc.add(acc,tr);
	  }
	  Res.set(i,j,acc);
	} 
     
      // Early zero detection 
      // --------------------
      // TODO : also put heuristic for global precision exhaustion 
      r1=0.0;
      for (j=0; j<d; j++) {
	tr.abs(Res.get(n-1,j)); // Change for several relations and ldim 
	if (tr.cmp(r1) > 0) r1=tr;
      }

      if (epsilon.cmp(r1) >0) {
       
	notfound=0; 
	// ICI 
	cout << endl << " WITH RESIDUE ENDING " << endl << endl; 
      }
 
      // Residue normalization  
      acc=0.0;
      for (i=0; i<n; i++)
	for (j=0; j<d; j++) {
	  tr.abs(Res(i,j));
	  if (tr.cmp(acc)>0) acc=tr;
	}
    
      for (j=0; j<d; j++)
	Res.div(j,0,Res.getcol(j,0),acc,n);

      for (j=0; j<d; j++)
	Resft.setcol(j,Res.getcol(j),0,n);

 
    } // End else not found, continue 

  } // End main loop on successive decompositions  

 

  // ICI 
  cout << endl << "   Part time: " << tps/1000 << " ms" << endl;
  cout << endl << "   Nb iterations " << nblov << endl;

  // The nullspace C should be in the last n+d-rel columns of V=Transpose(U^(-1))
  // See analogous procedure in nullspace indirect (on other entry types) 
  // ******************************************************************************
 
 

  ZZ_mat<mpz_t> UT;
  UT.resize(n+d,n+d);
  transpose(UT,V);
  
  C.resize(n,nbrel);
  
  for (i=0; i<n; i++)
    for (j=0; j<nbrel; j++) 
      C(i,j)=UT(n+d-nbrel+j,d+i); 
  
  
  // Check of the nullity 
  // Coul maybe not be used id decomp has succeeded?  
  // ***********************************************

  
  matrix<FP_NR<RT> > Cfp, ZM;

  Cfp.resize(n,nbrel);
  for (i=0; i<n; i++)
    for (j=0; j<nbrel; j++)
      set_z(Cfp(i,j),C(i,j));

  ZM.resize(d,nbrel);

  matprod(ZM,A,Cfp); // Should be epsilon-zero 
  
  /* Tuning of the epsilon = "zero" value for the test, about the floating point precision */
  /* ***************************************************************************************/

  
  FP_NR<RT> tmp1;

  /* Zero test    */
  /****************/
  
  int nullity=0;
  bool zeroedcol;

  for (j=0; j<nbrel; j++) {
       
    zeroedcol=true;

    for (i=0; i<d; i++) {
    
      tmp1.abs(ZM(i,j));     
      if (tmp1.cmp(epsilon) > 0)   
	zeroedcol=false;
    }
    if (zeroedcol) 
      nullity+=1;   
  }

 
  if (nullity != 1) {
    cout << " ***  ERROR: Input matrix non-full row rank or problem with the decomposition/precision" << endl; 
    cout << " ***          nullity = " << nullity << " <> " << nbrel << endl;  

    for (j=0; j<nbrel; j++) {
      
      tr=0.0;
      for (i=0; i<d; i++) {
	tmp1.abs(ZM(i,j));
	if (tmp1.cmp(tr) >0) tr=tmp1;
      }
    }
    cout << " ***          max residue = " << tr << "  >  epsilon = " << epsilon << endl;
    
    return -1;
  } 

  return(nbrel);
  
} // End body of FP relations HJLS  

} // end namespace hplll


#endif 
