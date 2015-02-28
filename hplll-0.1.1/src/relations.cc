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

  TO DO 

At least one >= 1 ? 

**************************************************************************************/ 

 
int relation_d(ZZ_mat<mpz_t>& C, const matrix<FP_NR<mpfr_t> >& A, long setprec) {

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
  found=relation_lift<mpz_t, double, matrix<FP_NR<double> > > (C, L, setprec, 0.99);
  
  
  start=utime()-start;

  
  cout << "   time internal: " << start/1000 << " ms" << endl;
  
  return found;

} 

  
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
