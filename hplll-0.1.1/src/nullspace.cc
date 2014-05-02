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



#ifndef HPLLL_NULLSPACE_CC
#define HPLLL_NULLSPACE_CC

namespace hplll { 

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

  F.decomp(1.1547005384, d); 
  cout << "nbswaps: " << F.nblov<< endl; 

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

  //print2maple(B,n,d); 

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
  Fgas<RT, mpz_t, RT,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, matrix<FP_NR<RT> > > G(B,NO_TRANSFORM,setprec);
  

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
  Fgas<mpfr_t, ZT, FT,  matrix<FP_NR<mpfr_t> >, MatrixZT, MatrixFT> F(FF, NO_TRANSFORM, setprec);
  

  // Decomposition 
  // -------------
  int flag_decomp;
  flag_decomp=F.decomp(1.1547005384, n-d);
  //flag_decomp=F.decomp(1.414213562, n-d);
  cout << "Nb Lovasz tests " << F.nblov << endl; 

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
  flag_decomp=F.decomp(1.1414, n-d);
  cout << "Nb tests " << F.nblov << endl; 

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

   int nullity=0;
   bool zeroedcol; 
   
   vector<int> nullcols;  // Zeroed columns in the upper part : 1, otherwise 0
   nullcols.resize(n); 
 

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

     if (lllmethod == L1) {

       ZZ_mat<mpz_t> tmpmat;
       tmpmat.resize(m,n);

       l1(tmpmat, B, 800, FPLLL);
      
       for (i=0; i<m; i++)
	 for (j=0; j<n; j++) 
	   B(i,j)=tmpmat(i,j); 


     } 
     if (lllmethod == HLLL) {
        
       Lattice<ZT, FT, MatrixZT, MatrixFT> D(B,NO_TRANSFORM,DEF_REDUCTION);
	 
       D.hlll(0.999999999973046945,LLL_VERBOSE);
       
       // Bad hack since no assignment for matrices ... 
       ZZ_mat<mpz_t> tmpmat;
       tmpmat.resize(n,m);
       transpose(tmpmat,D.getbase());
       transpose(B,tmpmat);

       //ICI 
       cout << endl << "    Nb swaps: " << D.nblov << endl ;
	 
     }
     if (lllmethod == FPLLL) {
	 
       ZZ_mat<mpz_t> AT;
       AT.resize(n,m);

       transpose(AT,B);

       //lllReduction(AT, 0.99, 0.62, LM_HEURISTIC,FT_DPE,0);
       lllReduction(AT, 0.75, 0.5, LM_WRAPPER,FT_DEFAULT,0);
       transpose(B,AT);
	
     }


     // Zero column test 
     // ----------------
     
     nullity=0;
  
     for (j=0; j<n; j++) {

       nullcols[j]=0; 
       zeroedcol=true;
       for (i=0; i<d; i++) if (B(i,j).sgn() !=0)  zeroedcol=false;  
       
       if (zeroedcol) { 
	 nullity+=1;
	 nullcols[j]=1; 
       } 
     }


     // Shift too small for being a unique shift 
     // ----------------------------------------
    
     if (nullity<n-d) { // attention dimensions 
       cout << " *** Anomaly, nullity is "  << nullity << " < " << n-d << ", probably a too small shift" << endl; 
       return -1; 
     }

   } // End unique shift 

   // Several shifts Lehmer like loop 
   // -------------------------------

   else {

     ZZ_mat<ZT> BL;
     BL.resize(m,n);

     bool notfound=1;
     
     Lattice<ZT, FT, MatrixZT, MatrixFT> Bt(B,TRANSFORM,DEF_REDUCTION);
     
     int max2up=0;
     max2up=maxbitsize(B,0,d,n);
     
     int current_shift;
     current_shift = -max2up;
 

     // Loop of elementary shifts 
     // -------------------------
     
    
     while (notfound) {


       BL=B;

       current_shift+=lsigma;
       // ICI        
       //cout << "current shift " << current_shift << endl; 
       
      shift_in(BL, current_shift, d);

      
      Bt.put(BL, d, lsigma+n); // Heuristic 
      

      Bt.hlll(0.999999999);
      
      // Required multiplication update when BL has been truncated 
      matprod_in(B,Bt.getU());
      // if not truncated one should not multiply evrything        

       
       // Test de colonnes à zéro 
       // -----------------------
       
       nullity=0;
   
       for (j=0; j<n; j++) {

	 nullcols[j]=0; 
	 zeroedcol=true;
	 for (i=0; i<d; i++) if (B(i,j).sgn() !=0)  zeroedcol=false;  
       
	 if (zeroedcol) { 
	   nullity+=1;
	   nullcols[j]=1; 
	 } 
       }

       //cout << "Nullity " << nullity << endl; 

       if (nullity>=n-d) {
	 notfound=0;
       } 
       
     } // End while loop on the Lehmer shifts 
     cout << "Shift " << current_shift << endl; 
   } // End else Lehmer case 



   // Output matrix 
   // -------------

   C.resize(n,nullity);

   int decj=0;
   for (j=0; j<n; j++) {  
     if (nullcols[j]==1) {
       for (i=0; i<n; i++) 
	 C(i,decj)=B(i+d,j);
       decj++; 
     }
   }


   return nullity;
};




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
       
     D.hlll(0.999999);
       
       B=D.getmixedbase();

       B.shift(-lalpha);

       

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


template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
nullspace_fp_hjls(ZZ_mat<ZT>& C, ZZ_mat<ZT> AZ, long int setprec) { 

  long d=AZ.getRows();
  long n=AZ.getCols();

  int nbrel;

  mpfr_set_default_prec(setprec);

  matrix<FP_NR<mpfr_t> > A;   // Input matrix 
  A.resize(d,n);
 
  int i,j;

  FP_NR<mpfr_t> tmp,maxa;

  maxa=0.0;
  for (i=0; i<d; i++) {
    for (j=0; j<n; j++) {
      set_z(tmp,AZ(i,j));
      tmp.abs(tmp);
      if (tmp.cmp(maxa) >0) maxa=tmp;
    }
  }


  for (i=0; i<d; i++) {
    for (j=0; j<n; j++) {
      set_z(tmp,AZ(i,j));
      tmp.div(tmp,maxa);
      A.set(i,j,tmp);
    }
  }

  //nbrel=relations_hjls<mpfr_t,dpe_t, MatrixPE<double, dpe_t> >(C,A,1,setprec);
  //nbrel=relations_hjls<mpfr_t,mpfr_t, matrix<FP_NR<mpfr_t> > >(C,A,n-d,setprec);
  nbrel=restarting_hjls<mpfr_t, double, matrix<FP_NR<double > > >(C,A,1,setprec);

  cout << "nbrel " << nbrel << endl;

  //print2maple(C,n,nbrel);
  return 0;
}	 

} // end namespace hplll


#endif 
