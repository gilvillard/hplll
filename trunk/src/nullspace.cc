/* Integer matrix nullspace  

Created Created Jeu  7 mar 2013 15:01:41 CET  
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


#include "decomp.cc"
#include "hlll.cc" 
#include "nullspace.h"

#ifndef HPLLL_NULLSPACE_CC
#define HPLLL_NULLSPACE_CC


// ********   ATTENTION   ****************************************
//
// Full rank, nullspace dimension a priori known 
//  Template RT but uses (may not be needed setprec of mpfr .....)
//
/* ***************************************************************/


// Faire une fct d'aiguillage générale avec un type method 

// Direct integer FGAS decomposition 
// *********************************


template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
nullspace_direct_decomp(ZZ_mat<ZT>& C, ZZ_mat<ZT> A) { 

  long d=A.getRows();
  long n=A.getCols();

  ZFgas<ZT, FT, MatrixZT, MatrixFT>F(A,TRANSFORM); 
 
  F.decomp(1.4142, d); 

  // Check of nullity and location of the zero columns 
  // ------------------------------------------------- 

  ZZ_mat<ZT> B;
  B.resize(n,d);
  transpose(B, F.getbase());

  int i,j,nullity=0;
  bool zeroedrow;
  vector<int> nullrows;
  nullrows.resize(n);   

  // Should be the first rows by construction in decomp 
  for (i=0; i<n; i++) {
       
    nullrows[i]=0; 
                      
    zeroedrow=true;
    for (j=0; j<d; j++) 
      if (B(i,j).sgn() !=0)  zeroedrow=false;  
       
    if (zeroedrow) { 
      nullity+=1;
      nullrows[i]=1; 
    } 
  }

  if (nullity != n-d) {
    cout << " *******  ERROR: Input matrix non-full row rank or problem with the decomposition" << endl; 
    cout << " *******          nullity = " << nullity << " <> " << n-d << endl;  
    return -1;
  } 

  // Output nullspace 
  // ****************
  
  C.resize(n,n-d);


  // Re-use of B 
  B.resize(n,n);
  transpose(B,F.getU());

  int decj=0;
  for (i=0; i<n; i++) {
    if (nullrows[i]==1) {
      for (j=0; j<n; j++) 
	C(j,decj)=B(i,j);
      decj++; 
    }
  }

  return(nullity);
 
} 



// FGAS generation + decomposition (PSLQ spirit)  
// *********************************************

// Not templated w.r.t. all matrices for the moment at this upper level
// Only the fgas generation 
 
// Note: const for A not possible since we use the Get function 
template<class RT, class FT, class MatrixFT> int  
fgasgen(matrix<FP_NR<RT> >& F,  ZZ_mat<mpz_t> A, const long setprec) { 


  long d=A.getRows();
  long n=A.getCols();

  int i,j;


  mpfr_set_default_prec(setprec); 

  // Fgas computation 
  // ----------------

  matrix<FP_NR<RT> > B;
  B.resize(n,n+d);

  FP_NR<RT> tmp;

  for (i=0; i<n; i++) {
    for (j=0; j<d; j++) {
      set_z(tmp,A(j,i));
      B.set(i,j,tmp);
    }
    B(i,d+i)=1;
  }

  F.resize(n-d,n);

  // Only for computing the orthogonalization, FT=RT 
  Fgas<RT, mpz_t, RT,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, matrix<FP_NR<RT> > > G(B,setprec);
  

  G.householder();

  // Re-use of B as a temp matrix 
  B=G.getR();

  // The fgas associated to A "à la" PSLQ 
  // ------------------------------------

  for (i=0; i<n-d; i++) 
    for (j=0; j<n; j++) 
      F(i,j)=B.get(d+i,d+j);
    
  return(0);
}


template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
nullspace_indirect_decomp(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, long int setprec=53) { 

  long d=A.getRows();
  long n=A.getCols();

  mpfr_set_default_prec(setprec); // May not be needed 

  matrix<FP_NR<mpfr_t> > FF;
  fgasgen<mpfr_t, FT, MatrixFT>(FF, A, setprec); 

  // FGAS
  // ----
  Fgas<mpfr_t, ZT, FT,  matrix<FP_NR<mpfr_t> >, MatrixZT, MatrixFT> F(FF,setprec);
  

  // Decomposition 
  // -------------
  int flag_decomp;
  flag_decomp=F.decomp(1.4142, n-d);

  if (flag_decomp < 0) return -1;

  // Check of nullity and location of the zero columns 
  // See analogous procedure in relations nullspace (on other entry types) 
  // and nullspace hjls below 
  // --------------------------------------------------------------------- 
  int i,j; 

  ZZ_mat<ZT> VT;
  VT.resize(n,n);
  transpose(VT,F.getV());

  C.resize(n,n-d);

  for (i=0; i<n; i++)
    for (j=0; j<n-d; j++) 
      C(i,j)=VT(d+j,i); 


  ZZ_mat<ZT> ZM(d,n-d);

  matprod(ZM,A,C);


  // Zero matrix test 
  
  int nullity=0;
  bool zeroedcol;

  for (j=0; j<n-d; j++) {
       
    zeroedcol=true;

    for (i=0; i<d; i++) 
      if ((ZM(i,j)).sgn() !=0) 
	zeroedcol=false;

    if (zeroedcol) 
      nullity+=1;   
  }

  
  if (nullity != n-d) {
    cout << " *******  ERROR: Input matrix non-full row rank or problem with the decomposition/precision" << endl; 
    cout << " *******          nullity = " << nullity << " <> " << n-d << endl;  
    return -1;
  } 


  return(nullity);

}



// Direct integer FGAS decomposition + HJLS dec 
// ********************************************


template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
nullspace_hjls(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, long int setprec=53) { 

  long d=A.getRows();
  long n=A.getCols();
 
  ZZ_mat<ZT> B;

  B.resize(n,n+d);

  int i,j;

  for (i=0; i<n; i++) {
    for (j=0; j<d; j++)
      B(i,j)=A(j,i);  // = is allowed, set is not needed? 
    B(i,d+i)=1;
  }
  
  ZFgas<ZT, FT, MatrixZT, MatrixFT>F(B,INV_T_TRANSFORM,d); 
  F.setprec(setprec);

  int flag_decomp;
  flag_decomp=F.decomp(1.4142, n-d);
  
  if (flag_decomp < 0) return -1;

  // Check of nullity and location of the zero columns 
  // See analogous procedure in relations nullspace (on other entry types) 
  // and nullspace hjls below 
  // --------------------------------------------------------------------- 
 

  ZZ_mat<ZT> VT;
  VT.resize(n+d,n+d);
  transpose(VT,F.getV());

  C.resize(n,n-d);

  for (i=0; i<n; i++)
    for (j=0; j<n-d; j++) 
      C(i,j)=VT(2*d+j,d+i); 


  ZZ_mat<ZT> ZM(d,n-d);

  matprod(ZM,A,C);


  // Zero matrix test 
  
  int nullity=0;
  bool zeroedcol;

  for (j=0; j<n-d; j++) {
       
    zeroedcol=true;

    for (i=0; i<d; i++) 
      if ((ZM(i,j)).sgn() !=0) 
	zeroedcol=false;

    if (zeroedcol) 
      nullity+=1;   
  }

  if (nullity != n-d) {
    cout << " *******  ERROR: Input matrix non-full row rank or problem with the decomposition/precision" << endl; 
    cout << " *******          nullity = " << nullity << " <> " << n-d << endl;  
    return -1;
  } 

  return(nullity);
}


/* ***********************************************

   Lehmer like LLL nullspace 

   d x n  input matrix 

   B m x n global lattice with the initial identity 

	  Returns the nullity as soon as >= n-d 
          Full rank required ? 

   alpha : bits, global shit for the upper part 
   lsigma : bits, elementary shift (step) size à la Lehmer 
   unique shift if lsigma=0;

   method : HLLL or FPLLL 

   ********************************************** */

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
nullspace_lll(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int lalpha, int lsigma=0, int lllmethod=HLLL) { 


  int d,n,m;
  int i=0,j=0;

   d=A.getRows();
   n=A.getCols();
   m=d+n;

   vector<int> nullcols;  // Zeroed columns in the upper part : 1, otherwise 0
   nullcols.resize(n); 
   for (j=0; j<n; j++) nullcols[j]=0;

   bool zeroedcol; 
   int nullity=0;  

   // Copy in the temporary basis and initial shift 
   // ---------------------------------------------

   ZZ_mat<ZT> B;

   B.resize(m,n);

   for (i=0; i<d; i++)
     for (j=0; j<n; j++) 
       B(i,j).mul_2exp(A(i,j),lalpha); 

   for (i=0; i<n; i++) 
       B(d+i,i)=1;  

 
   // Unique shift 
   // ------------
   if (lsigma ==0) {  

     if (lllmethod == HLLL) {
	
	 Lattice<ZT, FT, MatrixZT, MatrixFT> D(B,NO_TRANSFORM,DEF_REDUCTION);
	 
	 D.hlll(0.9);
       
	 // Bad hack since no assignment for matrices ... 
	 ZZ_mat<mpz_t> tmpmat;
	 tmpmat.resize(n,m);
	 transpose(tmpmat,D.getbase());
	 transpose(B,tmpmat);
	 
       }
       if (lllmethod == FPLLL) {
	 
	 ZZ_mat<mpz_t> AT;
	 AT.resize(n,m);

	 transpose(AT,B);

	 //lllReduction(AT, 0.99, 0.62, LM_HEURISTIC,FT_DPE,0,LLL_VERBOSE);
	 lllReduction(AT, 0.9, 0.62, LM_HEURISTIC,FT_DPE,0);

	 transpose(B,AT);
	 print2maple(B,m,n);
       }
     } // End unique shift 

     // ***** Lehmer like *****************
     /*
     else {

  //int max2up=0,max2down=0;
   
   // *********** Loop on the global shift for the nullspace computation
   //int K=0; 
   //int K2=0;
   
   //long w=0;
   //long maxw=1073741823;

       L1exp<mpz_t, double, dpe_t> Bt(B,true); //   <==== pas double 
   //L1<mpz_t, double> Bt(B,true);  // <== double 

  while (w < maxw) {

       // Bit lengths 
       // -----------
     
       // Of the m first rows  
       max2up=0;
       for (i=0; i<m; i++) // Attention aux dimensions 
	 for (j=0; j<d; j++)
	   max2up=max(max2up, (int) mpz_sizeinbase(B(i,j).getData(),2)); 
     
       //cout << endl << "max2up " << max2up << endl;

       // Of the last n-m rows 
       max2down=0; 
       for (i=m; i<n; i++)  // Attention aux dimensions 
	 for (j=0; j<d; j++)
	   max2down=max(max2down, (int) mpz_sizeinbase(B(i,j).getData(),2)); 
      
       //cout << endl << "max2down " << max2down << endl;
     
       // Multiplication by the shift or gradual lifting 
       // ----------------------------------------------

       // No explicit multiplication (lifting) 
       if ( (max2up-max2down) > lsigma) { 
	 
	 // Keep a lsigma shift in the upper part
	 // and the entire entries in the bottom part (size about max2down+lsigma) if with parameter 0 

	 //Bt.put(B,lsigma+2*d,max2up-max2down-lsigma);    //   <==== double 
	 Bt.put(B,0,max2up-max2down-lsigma);     // pas double 

       }
       // Continue the multiplication (shift) of the upper 
       // part for discovering the nullspace 
       else {
	 
	 // cout << " shift : " << lsigma-max2up+max2down << endl; 
	 
	 Bt.shift(B,m,lsigma-max2up+max2down); // m 
	 K2+=1;
       }
       
       //cout << endl << " ------------------------------------ " << endl << endl;

       //print2maple(Bt.getbase(),n,d);

       Bt.hlll(0.99);
       
       matprod(B,Bt.getU());
     } // end ok Lehmer not unique 
}
     */ 

     // Test de colonnes à zéro 
     // -----------------------
     /*
     // Attention aux dimensions 
     nullity=0;

     for (j=0; j<d; j++) {

       nullcols[j]=0; // Put 0 at each time, there can be a permutation after the zero col is found
                      // A zero column may change (close to the end) 
       zeroedcol=true;
       for (i=0; i<m; i++) if (B(i,j).sgn() !=0)  zeroedcol=false;  
       
       if (zeroedcol) { 
	 nullity+=1;
	 nullcols[j]=1; 
       } 
     }

     if (nullity>=d-m) {
       K=w+1; // Kept for output 
       w=maxw+1;
     } 

     // Shift to small for being a unique shift 
     if (unique ==1) {
       if (nullity<d-m) { // attention dimensions 
	 cout << endl << " **********    Anomaly, nullity is "  << nullity << " < " << d-m << "  probably a too small shift" << endl; 
	 return -1; 
       }
     }

     w++;  
   } // Main loop on K, on the shifts for finding the nullspace 

   if (method ==1) 
     cout << endl << " ---- HLLL -----  Lovasz's tests: " << nbtests << "    Swaps: " << swaps << "  ------------ " << endl  <<  endl; 
   if (method ==2) 
     cout << endl << " ---- FPLLL -------------------------------- " << endl << endl; 

   if (unique==1) {
     cout <<  "Unique explicit initial shift of " << alpha << " bits" << endl ; 
   }
   if (unique ==0) {
     cout <<  K << " implicit or explicit shifts of " << lsigma << " bits each" << "  (alpha = " << alpha << ")" << endl ; 
     cout <<  K2 << "  explicit ones in the second phase" << endl ; 
   }

   // Output matrix 
   C.resize(d,nullity);

   int decj=0;
   for (j=0; j<d; j++) {  // Attention dimensiosn 
     if (nullcols[j]==1) {
       for (i=m; i<n; i++) 
	 C(i-m,decj)=B(i,j);
       decj++; 
     }
   }
     */



   // Stats 
   //******
   /*
   ZZ_mat<mpz_t> U;

   U.resize(m,d); 

   for (j=0; j<d; j++)  // Attention dimensions 
     for (i=0; i<m; i++) 
	 U(i,j)=B(i,j);
      
   cout << "Size upper U residue " << maxbitsize(U)-lalpha << endl;

   U.resize(d,d);
   for (j=0; j<d; j++) 
     for (i=0; i<d; i++) 
       U(i,j)=B(i+m,j);

   cout << "Size of unimodular U " << maxbitsize(U) << endl;
   */

   return nullity;
};


#endif 
