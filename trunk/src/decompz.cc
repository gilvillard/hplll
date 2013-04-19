/* FGAS decomposition 

Created Mar 18 jan 2011 18:10:25 CET 
Main update Jeu  7 mar 2013 14:14:54 CET
Copyright (C) 2011, 2012, 2013      Gilles Villard 

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


#ifndef HPLLL_DECOMPZ_CC
#define HPLLL_DECOMPZ_CC

#include "decompz.h"

// ********   ATTENTION   **********************

// Should be full rank on the left 
//
// Todo : proprifications, which ones and R update ? 
//
//=========================================================================
//
//  Global variable dec used for switching between HJLS and PSLQ
//   dec is the dimension e.g. dec=1 one vector of the space 
//   for simultaneous relations detection 
//   Which use in the integre case ? 
//  
//  (Assume that the first dec columns of the fgas are, d x dec 
//  are for input columns for searching relations 
//  The rest of the fgas is the identity (or something else) for HJLS)
// 
//
// FGAS Decomposition following the paper 
//
// Here for a mpfr matrix  
// Precision implicitly given from the construction hence the input FGAS 
//
// Should be extended for stopping at a given solution size 
// Target dim : dimension of the lattice part 
// Hence in the nullspace case should be lower than the input rank 
//
// The transformation applied to F is U where in the direct case 
// one finds the nullspace 
//
//==========================================================================
 
// ******  Difference with decomp.cc : computations ans tests on the residue Y 
//         and the relation size 
// ******     Computation of U instead of V for the indirect real case 

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
ZFgas<ZT, FT, MatrixZT, MatrixFT>::decomp(long double gamma, long int targetdim) { 

  // Verbose 
  //cout << "Direct mpfr (floating point) FGAS decomposition ...."  << endl;

  int i,j;

  vector<FP_NR<FT> >  g(d);
  g[0]=gamma;
  for (i=1; i<d; i++) g[i].mul(g[i-1],g[0]); 

  FP_NR<FT> tmp1,tmp2; 
  
  int maxh=dec;

  FP_NR<FT> diagold; // For the gap test 

  int prevkappa=-1; // For the looping test betwenn tow indices 
  vector<FP_NR<FT> >  prevR(d);
 

  // Initial size reduction and implicit Householder 
  // ***********************************************

  for (j=0; j<n; j++) hsizereduce(j);   // Implicitly computes Householder 


  // Initialization for the confidence gap value
  // minimum of the diagonal R 
  // ******************************************* 
  
  diagold.abs(R.get(dec,dec));
  for (i=dec+1; i<d; i++) {
    tmp1.abs(R.get(i,i));
    if (tmp1 < diagold)
      diagold=tmp1;
  }
  

  // **********************************************************
  // Main loop 
  // The update of R is done at the end of the work in the loop
  // **********************************************************
  
  int start;

  start=utime();
   
 
  while (true) {

    nblov++;

    if (((nblov%800000)==0) && (nblov > 0))   cout << nblov << " tests" << endl;  

    // For the gap test 
    // ****************
    
    tmp1.div(R.get(d-ldim-1,d-ldim-1),diagold);
    tmp1.abs(tmp1);

    // tmp1 kept for the nextcoming test 

    // New dimension test for the lattice component 
    // ******************************************** 
 
    // Heuristic test,  purely input data dependent (no link with the FT précision?)  
    // Anyway, problem here detected later: the "relation" will not be a relation 
    if (tmp1.cmp(confidence_gap) < 0) { 

      
      // We put the non-zero column at the end for respecting the decomp structure 
      if (maxh < n-ldim-1) {
    
	F.colswap(maxh,n-ldim-1);

	if (transf)
	  U.colswap(maxh,n-ldim-1);
	else 
	  W.colswap(maxh,n-ldim-1);

	// No need of swapping in R that will be updated later  
	// Note that the swap here is in the right rectangular part 
      }

      ldim+=1;   
           
    }

    // Global termination test and stopping heuristics
    // ***********************************************

    // Can be put in the while condition
    if (ldim >= targetdim) {  // Termination could be anticipated a few by testing zero results 
      
      return(0);
    }

    // For the gap test, check in R 
    // ----------------------------
   
    // In case dec = 0, dec could be put to 0 
    diagold.abs(R.get(dec,dec));  // Test only R.get(d-ldim-1,d-ldim-1) for efficiency?
    for (i=dec+1; i<d-ldim-1; i++) {
      tmp1.abs(R.get(i,i));
      if (tmp1 < diagold)
	diagold=tmp1;
    }
  

    // Swap position search in the left square part 
    // ******************************************** 

    maxh=dec;
    tmp1.mul(g[0],R.get(dec,dec));
    tmp1.abs(tmp1);

    for (j=dec+1; j<d-ldim; j++) {
   
      tmp2.mul(g[j],R.get(j,j));
      tmp2.abs(tmp2);
      
      if (tmp2 > tmp1) {
	maxh=j;    // (si plusieurs identiques ?)
	tmp1=tmp2;
      }
    }  


    // Heuristic check: precision exhausted, some invariant not preserved 
    // See analogous check in hlll
    // ******************************************************************

    if (maxh<d-ldim) {

      if(prevkappa==maxh) {

	FP_NR<FT> t;
	t.abs(R.get(maxh,maxh));
        if (t > prevR[maxh]) { 

	  cout << " **** #tests = " << nblov << "  Anomaly, norm increased, the FT precision may be exhausted" << endl; 
	  return -1;
	}
      }

      prevkappa=maxh;
      prevR[maxh].abs(R.get(maxh,maxh));
    }

    // Swap ok in the left part 
    // ************************  

    if (maxh<d-ldim-1) {
      
      
      F.colswap(maxh,maxh+1);

      if (transf)      
	U.colswap(maxh,maxh+1);
      else 
	W.colswap(maxh,maxh+1);

      // Restricted to the left part, right part below in the else   
      for (j=maxh; j<d-ldim; j++) hsizereduce(j);   // Householder implicitly computed 

    }

    // Swap position search in the right rectangular part 
    // ************************************************** 

    else { 

      // On ne fait Householder et la proprification que si on 
      // rentre dans la zone de recherche 

      for (j=d-ldim; j<n-ldim; j++) hsizereduce(j);   // Householder implicitly computed 

      maxh=d-ldim;
      tmp1.abs(R.get(d-ldim-1,d-ldim));

      for (j=d-ldim+1; j<n-ldim; j++) {
	
	tmp2.abs(R.get(d-1-ldim,j));
 
	if (tmp2 > tmp1) {
	  maxh=j;
	  tmp1 = tmp2;
	}

      } 

      F.colswap(d-ldim-1,maxh);

      if (transf)
	U.colswap(d-ldim-1,maxh);
      else 
	W.colswap(d-ldim-1,maxh);

      // Needed: this column can be used for swap if non-zero diag at next phase 
      // Not reduced everywhere for not reducing w.r.t. a zero (very small) entry 
      //  the rest will be reduced at next phase in the right part if non-zero diag 

      hsizereduce(d-ldim-1);

    }  // End else right part 
    
  } // End main iteration loop 

 return 0; 

} 


/* -------------------------------------------------------------------------
   Size reduction 

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm
   ------------------------------------------------------------------------- */
// exactly same as in decomp.cc but FP_NR<RT> xr; not necessary here 
// hence F with xz 

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
ZFgas<ZT, FT, MatrixZT, MatrixFT>::hsizereduce(int kappa) { 


  FP_NR<FT> approx;
  
  approx=0.01;

  FP_NR<FT> x,t,tmpfp;
  Z_NR<ZT>  xz,tmpz;

  int i,w=0;

  bool nonstop=1;
  bool somedone=0;

  int elim_index;
   
  householder(kappa); 
 
  while (nonstop) {

    w++;

    somedone = 0;

    elim_index=min(kappa-1,d-1-ldim);

    
    for (i=elim_index; i>-1; i--){
   
      x.div(R.get(i,kappa),R.get(i,i)); 

      x.rnd(x);

      if (x.sgn() !=0) { 

	set_f(xz,x);
	
	if (xz == -1){
	  somedone=1;

	  R.addcol(kappa,i,i);

	  F.addcol(kappa,i,d);

	  if (transf)
	    U.addcol(kappa,i,n);
	  else 
	    W.subcol(i,kappa,n);

	}
	// ----------------------------------------------- 
	else if (xz == 1){
	  somedone=1;

	  R.subcol(kappa,i,i);

	  F.subcol(kappa,i,d);
	  if (transf)
	    U.subcol(kappa,i,n);
	  else
	    W.addcol(i,kappa,n);
	}
	// ------------------------------------------------ 
	else {   
	  somedone=1;

	  R.submulcol(kappa,i,x,i);
	  F.submulcol(kappa,i,xz,d);

	  //if (DPE_EXP(x.getData()) < 31) {

	    //int ii;
	    //ii=xz.get_si();

	    //if (ii > 0) 
	      //vector_addmul_ui(U.getcol(i),U.getcol(kappa),ii,n);
	    //else 
	      //vector_submul_ui(U.getcol(i),U.getcol(kappa),-ii,n); 
	  //}

	  //else {
	  if (transf)
	    U.submulcol(kappa,i,xz,n);
	  else 
	    W.addmulcol(i,kappa,xz,n);
	    
	}
	// ----------------------------------------------- 

      } // End non zero combination 

    } // End loop through the column
    

    if (somedone) {

      t.mul(approx,normB2[kappa]);

      householder(kappa);  // Implicit update of normB2 
    
      nonstop = (normB2[kappa] < t);  // Some decrease ? 

    }
    else {
      nonstop=0;
    }
  } // end while 

  return 0;

};


/* --------------------------------------------- */
/* Householder on column kappa                   */
/*     partially : d-ldim x  n                   */
/* --------------------------------------------- */
// exactly same as in decomp.cc

// Change w.r.t. hlll _up and _v and we consider ldim 

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
ZFgas<ZT, FT, MatrixZT, MatrixFT>::householder(int kappa) {

  int k;
  FP_NR<FT> nrtmp,s,w; 
    
  // Left non-singular part
  // ----------------------
  if (kappa < d-ldim)  {  // Introduction of ldim 
    
    
    R.setcol(kappa,F.getcol(kappa),0,d); 
    
    fp_norm_sq(normB2[kappa], R.getcol(kappa), d);

    for (k=0; k<kappa; k++) {
      
      scalarprod(nrtmp, V.getcol(k,k), R.getcol(kappa,k), d-k);
      
      R.fmasub(kappa,k,R.getcol(kappa,k), V.getcol(k,k), nrtmp, d-k); 
      
    }

    w=R.get(kappa,kappa);
    
    if (w >=0) {
      fp_norm(s,R.getcol(kappa,kappa),d-kappa);
      nrtmp.neg(s);
      R.set(kappa,kappa,nrtmp);  
    }
    else {
      fp_norm(nrtmp,R.getcol(kappa,kappa),d-kappa); 
      R.set(kappa,kappa,nrtmp);
      s.neg(nrtmp); 
    }
    
    w.add(w,s);
    s.mul(s,w);
    s.sqrt(s);
    
    V.div(kappa,kappa+1, R.getcol(kappa,kappa+1), s, d-kappa-1);
    nrtmp.div(w,s);
    V.set(kappa,kappa,nrtmp); 
        
    for(int i=kappa+1; i<d; i++) R.set(i,kappa,0.0); 
       
  }  // end of kappa in left part  

  else {

    // Non square part 
    // ---------------
    
     
    R.setcol(kappa,F.getcol(kappa),0,d);
    
    fp_norm_sq(normB2[kappa], R.getcol(kappa), d);
    
    for (k=0; k<d-ldim; k++) {  // Introduction of ldim 
      
      scalarprod(nrtmp, V.getcol(k,k), R.getcol(kappa,k), d-k);
	
      R.fmasub(kappa,k,R.getcol(kappa,k), V.getcol(k,k), nrtmp, d-k); 
      
    }   
  } // end of kappa in right part 

  return 0; 
}


/* --------------------------------------------- */
/* Complete Householder  d x  n                  */
/* --------------------------------------------- */
// exactly same as in decomp.cc

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
ZFgas<ZT, FT, MatrixZT, MatrixFT>::householder() {

  int k,kappa;
  FP_NR<FT> nrtmp,s,w; 
      
  // Full rank left square part 
  // --------------------------

  for (kappa=0; kappa<d; kappa++) {  

    R.setcol(kappa,F.getcol(kappa),0,d);

    for (k=0; k<kappa; k++) {
      
      scalarprod(nrtmp, V.getcol(k,k), R.getcol(kappa,k), d-k);      
      R.fmasub(kappa,k,R.getcol(kappa,k), V.getcol(k,k), nrtmp, d-k); 
	
    }

    w=R.get(kappa,kappa);
  
    if (w >=0) {
      fp_norm(s,R.getcol(kappa,kappa),d-kappa);
      nrtmp.neg(s);
      R.set(kappa,kappa,nrtmp);  
      }
    else {
      fp_norm(nrtmp,R.getcol(kappa,kappa),d-kappa);
      R.set(kappa,kappa,nrtmp);
  
      s.neg(nrtmp); 
    }
      
    w.add(w,s);
    s.mul(s,w);
    s.sqrt(s);

    V.div(kappa,kappa+1, R.getcol(kappa,kappa+1), s, d-kappa-1);
    nrtmp.div(w,s);
    V.set(kappa,kappa,nrtmp); 
        
    for(int i=kappa+1; i<d; i++) R.set(i,kappa,0.0); 
       
  }  // end left part 

  // Non square part 
  // ---------------

  for (kappa=d; kappa<n; kappa++) {

    R.setcol(kappa,F.getcol(kappa),0,d);

    for (k=0; k<d; k++) {

      scalarprod(nrtmp, V.getcol(k,k), R.getcol(kappa,k), d-k);
      R.fmasub(kappa,k,R.getcol(kappa,k), V.getcol(k,k), nrtmp, d-k); 

    }
  } // end right part 


  return 0; 
}


// ********************************
// VARIOUS   ACCESS CONSTRUCTION  
// ********************************

/*
 template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> 
   inline unsigned int ZFgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::getprec_fgas() { 

   // Actually the global ambiant mpfr prec (see also setprec_internal) 
   // FP_NR getprec is global, not the one of the element with mpfr 
  return (F.get(0,0)).getprec(); 
  
}
*/


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline unsigned int 
ZFgas<ZT,FT, MatrixZT, MatrixFT>::setprec(unsigned int prec) {

  // Re-initializations
  // Should be for each variable, non global? 
  // push and pop ? 
  // Put clear in the resize ? 

  unsigned oldprec;
  oldprec=getprec();

  mpfr_set_default_prec(prec);

  R.clear();
  R.resize(d,n);

  unsigned newprec;
  newprec=(R.get(0,0)).getprec(); 
  if (newprec == oldprec) cout << "Warning: in function setprec decompz, the change of precision has no effect" << endl; 

  V.clear();
  V.resize(d,d);

  normB2.clear();
  normB2.resize(n); 

  // The old precision 
  return oldprec;
}



template<class ZT,class FT, class MatrixZT, class MatrixFT> inline unsigned int 
ZFgas<ZT,FT, MatrixZT, MatrixFT>::getprec() {

  return (R.get(0,0)).getprec(); 
  
}

template<class ZT, class FT, class MatrixZT, class MatrixFT> 
inline  int ZFgas<ZT, FT, MatrixZT, MatrixFT>::set_confidence_gap(double ratio) {

  confidence_gap=ratio;

  return 1;
}


template<class ZT, class FT, class MatrixZT, class MatrixFT> 
inline ZZ_mat<ZT> ZFgas<ZT, FT, MatrixZT, MatrixFT>::getU() { 

  ZZ_mat<ZT> UU;
  
  if (transf) {

    UU.resize(n,n);
    for (int i=0; i<n; i++) 
      for (int j=0; j<n; j++) UU.Set(i,j,U(i,j)); // reprendre boucle sur les colonnes 
    
    return UU;
  }
  else {

    cout << "*** Access error: U has not been computed" << endl;
    return UU;

  }
}

template<class ZT, class FT, class MatrixZT, class MatrixFT> 
inline ZZ_mat<ZT> ZFgas<ZT, FT, MatrixZT, MatrixFT>::getV() { 

  ZZ_mat<ZT> WW;

  if (transf) {

    cout << "*** Access error: V has not been computed" << endl;
    return WW;

  }
  else {

    WW.resize(n,n);
    for (int i=0; i<n; i++) 
      for (int j=0; j<n; j++) WW.Set(i,j,W(i,j)); // reprendre boucle sur les colonnes 
    
    return WW;
  }
}



template<class ZT, class FT, class MatrixZT, class MatrixFT> 
inline ZZ_mat<ZT> ZFgas<ZT, FT, MatrixZT, MatrixFT>::getbase() { 

  ZZ_mat<ZT> B(d,n);
  for (int i=0; i<d; i++) 
    for (int j=0; j<n; j++) B.Set(i,j,F(i,j)); // reprendre boucle sur les colonnes 

  return B;
}



// Initialization and construction 
// -------------------------------

template<class ZT, class FT, class MatrixZT, class MatrixFT> void  
ZFgas<ZT, FT, MatrixZT, MatrixFT>::init(int d, int n, int forUV, long int inputdec) {


  if (forUV == TRANSFORM) transf=1;  // U 
  else transf=0;  // V 

  dec=inputdec;

  tmpcompt=0;

  nblov=0;

  ldim=0;

  // Inverse transpose of the transformation matrix (hence possibly nullspace) 
  if (transf) {  
    U.resize(n,n);
    for (int i=0; i<n; i++) U(i,i)=1; 
  }
  else {
    W.resize(n,n);
    for (int i=0; i<n; i++) W(i,i)=1; 
  }
  // !!! Put back the old precision if needed outside 
  // Pas forcément MPFR !!! 
  //  mpfr_set_default_prec(setprec);

  // Verbose 
  //cout << endl << "Decomp precision for mpfr: " <<  mpfr_get_default_prec() << " bits" << endl ;

  F.resize(d,n);
 
  R.resize(d,n);

  V.resize(d,d);

  normB2.resize(n);


  // Ratio non zero / zero 
  confidence_gap = 1e-9; 
 
}

// Construction 
template<class ZT, class FT, class MatrixZT, class MatrixFT>  
ZFgas<ZT, FT, MatrixZT, MatrixFT>::ZFgas(ZZ_mat<ZT> Finput, int forUV, long int inputdec=0) {

  d=Finput.getRows();
  n=Finput.getCols();

  //n=d-m;

  init(d,n, forUV, inputdec); 

  int i,j;

  for (i=0; i<d; i++) 
    for (j=0; j<n; j++) {
      F.set(i,j,Finput(i,j));
      } 

}


#endif
