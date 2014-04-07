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


#ifndef HPLLL_DECOMP_CC
#define HPLLL_DECOMP_CC

#include "decomp.h"

// ********   ATTENTION   **********************
//
// Nullity test ? Especially with long double ? 
//
// Should be full rank on the left 
//
// Attention !!! Change of precision, put back the previous one if 
// needed since this is done through  mpfr_set_default_prec(setprec); 
//
// Todo : proprifications, which ones and R update ? 
//
//=========================================================================
//
//  Parameter dec used for switching between HJLS and PSLQ
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
// The transformation applied to F is V: we compute U = Transpose (V^(-1)) 
//
//==========================================================================
 

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> int  
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::decomp(long double gamma, long int targetdim, long int dec=0) { 

  cout << "Direct mpfr (floating point) FGAS decomposition ...."  << endl;

  int i,j;

  vector<FP_NR<FT> >  g(d);
  g[0]=gamma;
  for (i=1; i<d; i++) g[i].mul(g[i-1],g[0]); 

  FP_NR<FT> tmp1,tmp2; // For determining the swap 

  int maxh=dec;

  FP_NR<FT> diagold; // For the epsilon-test 
  

  // Initial size reduction and implicit Householder 
  // ***********************************************
  
  int redmin; // Swap index for subsequent size reductions

  for (j=0; j<n; j++) hsizereduce(j);   // Implicitly computes Householder 


  diagold=R.get(dec,dec);

  // Initialization for the epsilon test value 
  for (i=dec+1; i<d; i++) {
    tmp1.abs(R.get(i,i));
    if (tmp1 > diagold)
      diagold=tmp1;
  }
  

  // **********************************************************
  // Main loop 
  // The update of R is done at the end of the work in the loop
  // **********************************************************
  
  int start;

  start=utime();
   
  while (true) {

    // ICI 
    if (nblov==2000) print2maple(Y,d,1); 
    if (nblov==3000) print2maple(Y,d,1); 
    if (nblov==4000) print2maple(Y,d,1); 

    // ICI 
    if (nblov==3470) {

      for (j=0; j<n; j++)
	Wf.setcol(j,W.getcol(j),0,n); 

      MatrixRT testmat;
      testmat.resize(n,n);
      matprod(testmat,F0,Wf);


      FP_NR<RT> tmp;

      for (i=0; i<d; i++) 
	for (j=0; j<n; j++) {
	  tmp.sub(testmat.get(i,j),F.get(i,j));
	  tmp.abs(tmp);
	  testmat.set(i,j,tmp);
	}

      FP_NR<RT> minf;
      FP_NR<RT> mrel;
      minf=0.0;
      mrel=0.0;

      for (i=0; i<d; i++) 
	for (j=0; j<n; j++) {
	  if (((testmat.get(i,j)).sgn() !=0) && ((Wf.get(i,j)).sgn() !=0)) {
	    
	    tmp.abs(testmat.get(i,j));
	    if (tmp.cmp(minf) > 0) minf=tmp;
	    tmp.div(testmat.get(i,j),Wf.get(i,j));
	    if (tmp.cmp(mrel) > 0) mrel=tmp;
	  }
	}
      cout << "****** Absolute error: ";
      minf.print();
      cout << "   Relative error: ";
      mrel.print();
      cout << endl; 

    } // End size test 

    nblov++;

    if (((nblov%800000)==0) && (nblov > 0))   cout << nblov << " tests" << endl;  

    // For the epsilon test 
    // ********************

    tmp1.div(R.get(d-ldim-1,d-ldim-1),diagold);
    tmp1.abs(tmp1);

    // New dimension test for the lattice component 
    // ******************************************** 

    //cout << "++++++ " << tmp1 << "       "  <<  R.get(d-ldim-1,d-ldim-1) <<  endl; 
 
    if (tmp1 < 1e-9) {  // ***** Voir comment faire ce TEST ? Purely algorithmic, no link with the FT précision ?  
                        // Anyway, problem here detected later: the "relation" will not be a relation 
     // Verbose 
     //cout <<  " ********  LDIM: " << ldim +1 <<  " maxh+1: " << maxh+1 << "    n-ldim: "  << n-ldim <<  "    " << tmp1 << endl; 
     
      // We put the non-zero column at the end for respecting the decomp structure 
      if (maxh < n-ldim-1) {
    
	F.colswap(maxh,n-ldim-1);
	U.colswap(maxh,n-ldim-1);
	W.colswap(maxh,n-ldim-1);
	{
	  FP_NR<RT> ty;
	  ty=Y.get(n-ldim-1-1,0);
	  Y.set(n-ldim-1-1,0,Y.get(maxh-1,0));
	  Y.set(maxh-1,0,ty);
	}
	//Y.rowswap
	// Pas besoin de swap pour R qui est mise à jour progressivement ensuite 
	// Note that the swap in the right rectangular part 
      }

      ldim+=1;         
    }


    // Global termination test 
    // ***********************
   // Can be put in the while condition 

    if (ldim >= targetdim) {  // Termination could be anticipated a few by testing zero results 
      // ICI       
      print2maple(Y,d,1);

      return(0);
    }

    // for the epsilon test 
    diagold=R.get(d-ldim-1,d-ldim-1);
    diagold.abs(diagold);


    // ICI 
    {
      FP_NR<FT> tt1,tt2;
      tt1.abs(R.get(dec,dec));
      for (j=dec+1; j<d-ldim; j++) {
	tt2.abs(R.get(j,j));
	if (tt2.cmp(tt1) > 0) tt1=tt2;
      }
      tt2=1.0;
      tt1.div(tt2,tt1); 

      Z_NR<ZT> res;
      set_f(res,tt1);
      cout << "*******  Relation > ";
      res.print();
      cout << endl; 


      ZZ_mat<mpz_t> C;
      C.resize(n,n);
      for (i=0; i<n; i++) 
	for (j=0; j<n; j++) 
	  C(i,j)=W.get(i,j);

      int size;
      size=maxbitsize(C);
      cout << "         W-size = " << size << endl; 

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



    // Swap ok in the left part 
    // ************************  

    if (maxh<d-ldim-1) {

      F.colswap(maxh,maxh+1);
      U.colswap(maxh,maxh+1);
      W.colswap(maxh,maxh+1);
      {
	FP_NR<RT> ty;
	ty=Y.get(maxh,0);
	Y.set(maxh,0,Y.get(maxh-1,0));
	Y.set(maxh-1,0,ty);
      }

      redmin=maxh;

      // On se limite à d-ldim, ci-dessous si passage dans l'autre zone  
      for (j=redmin; j<d-ldim; j++) hsizereduce(j);   // Effectue aussi Householder

    }

    // Swap position search in the right rectangular part 
    // ************************************************** 

    else { 

      // On ne fait Householder et la proprification que si on 
      // rentre dans la zone de recherche 

      for (j=d-ldim; j<n-ldim; j++) hsizereduce(j);   // Effectue aussi Householder 

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
      U.colswap(d-ldim-1,maxh);
      W.colswap(d-ldim-1,maxh);
      {
	FP_NR<RT> ty;
	ty=Y.get(d-ldim-1-1,0);
	Y.set(d-ldim-1-1,0,Y.get(maxh-1,0));
	Y.set(maxh-1,0,ty);
      }

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


template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> int  
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::hsizereduce(int kappa) { 


  FP_NR<FT> approx;
  
  approx=0.01;

  FP_NR<FT> x,t,tmpfp;
  Z_NR<ZT>  xz,tmpz;

  FP_NR<RT> xr;

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
	  U.subcol(i,kappa,n);
	  W.addcol(kappa,i,n);
	  if (i>=1) {
	    FP_NR<RT> ty;
	    ty=0.0;
	    ty.sub(Y.get(i-1,0),Y.get(kappa-1,0));
	    Y.set(i-1,0,ty);
	  }
	}
	/* ----------------------------------------------- */
	else if (xz == 1){
	  somedone=1;

	  R.subcol(kappa,i,i);

	  F.subcol(kappa,i,d);
	  U.addcol(i,kappa,n);
	  W.subcol(kappa,i,n);
	  if (i>=1) {
	    FP_NR<RT> ty;
	    ty=0.0;
	    ty.add(Y.get(i-1,0),Y.get(kappa-1,0));
	    Y.set(i-1,0,ty);
	  }
	  
	}
	/* ----------------------------------------------- */
	else {   
	  somedone=1;

	  set_z(xr,xz);

	  R.submulcol(kappa,i,x,i);
	  F.submulcol(kappa,i,xr,d);

	  /*if (DPE_EXP(x.getData()) < 31) {

	    int ii;
	    ii=xz.get_si();

	    if (ii > 0) 
	      vector_addmul_ui(U.getcol(i),U.getcol(kappa),ii,n);
	    else 
	      vector_submul_ui(U.getcol(i),U.getcol(kappa),-ii,n); 
	  }

	  else {*/
	  
	  vector_addmul(U.getcol(i),U.getcol(kappa),xz,n);
	  W.submulcol(kappa,i,xz,n); 

	  if (i>=1) {
	    FP_NR<RT> ty;
	    ty=0.0;
	    ty.mul(xr,Y.get(kappa-1,0));
	    ty.add(Y.get(i-1,0),ty);
	    Y.set(i-1,0,ty);
	  }

	}
	/* ----------------------------------------------- */

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

// Change w.r.t. hlll _up and _v and we consider ldim 

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> int  
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::householder(int kappa) {

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

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> int  
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::householder() {

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


template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> 
inline ZZ_mat<ZT> Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::getU() { 

  ZZ_mat<ZT> UU(n,n);
  for (int i=0; i<n; i++) 
    for (int j=0; j<n; j++) UU.Set(i,j,U(i,j)); // reprendre boucle sur les colonnes 

  return UU;
}


// Initialization and construction 
// -------------------------------

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> void  
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::init(int d, int n, long setprec) {

  tmpcompt=0;

  nblov=0;

  ldim=0;

  // Inverse transpose of the transformation matrix (hence possibly nullspace) 
  U.resize(n,n);
  for (int i=0; i<n; i++) U(i,i)=1; 

  // Transformation matrix 
  W.resize(n,n);
  for (int i=0; i<n; i++) W(i,i)=1; 

  // !!! Put back the old precision if needed outside 
  mpfr_set_default_prec(setprec);

  // Verbose 
  //cout << endl << "Decomp precision for mpfr: " <<  mpfr_get_default_prec() << " bits" << endl ;

  F.resize(d,n);
  F0.resize(d,n);
  Wf.resize(n,n);

  Y.resize(d,1);

  R.resize(d,n);

  V.resize(d,d);

  normB2.resize(n);
}

// Construction 
template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>  
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::Fgas(MatrixRT Finput, long setprec) {

  d=Finput.getRows();
  n=Finput.getCols();

  //n=d-m;

  init(d,n,setprec); 

  int i,j;

  for (i=0; i<d; i++) 
    for (j=0; j<n; j++) {
      F(i,j)=Finput(i,j);
      } 

  for (i=0; i<d; i++) 
    for (j=0; j<n; j++) {
      F0(i,j)=Finput(i,j);
    } 

  for (i=0; i<d; i++) 
    Y(i,0)=Finput(i,0);

}


#endif
