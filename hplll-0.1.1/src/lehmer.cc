/* Integer matrix nullspace  

Created Created Sam 11 mai 2013 15:52:10 CEST   
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
#include "ratio.h" 

#ifndef HPLLL_LEHMER_CC
#define HPLLL_LEHMER_CC

// METTRE DANS LE NAMESPACE 
using namespace hplll; 



/* ***********************************************

   Lehmer like LLL  

   (d+m) x n  input matrix 

   ********************************************** */


  
template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
lehmer_f(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int lsigma, double delta) { 

  if (lsigma==0) {
    
    Lattice<ZT, FT, MatrixZT, MatrixFT> B(A,NO_TRANSFORM,DEF_REDUCTION); 
    
    B.hlll(delta);
    
    C=B.getbase();
    
  }
  // Non trivial shift 
  // -----------------
  
  else {   // lsigma <> 0 

    int i,j;
    
    int m=1;   // To change in parameter

    int d=A.getCols();

    C.resize(m+d,d);

    for (i=0; i<m+d; i++)
       for (j=0; j<d; j++)
	 C(i,j)=A(i,j);
    
    int bitsize = maxbitsize(A,0,m,d);

    int def = -bitsize;

    ZZ_mat<double> Af;
    Af.resize(m+d,d);

    Z_NR<ZT> tz;

    Lattice<double, double,  matrix<Z_NR<double> >, matrix<FP_NR<double> > > B(Af,TRANSFORM,DEF_REDUCTION);

    ZZ_mat<long> U;
    U.resize(d,d);
    
    ZZ_mat<double> Uf;
    Uf.resize(d,d);

    FP_NR<double> tf;
     
    // Main Lehmer loop on the defect
    // ------------------------------
    
    while (def < 0) { 

      def = min(0,def+lsigma);

      // ICI 
      int mmax=0;
      for (j=0; j<d; j++)
	if (mmax <  size_in_bits(C(0,j))) mmax=size_in_bits(C(0,j));
	
      def=lsigma-mmax;

      cout << endl << "def: " << def << endl; 

	
      
      for (i=0; i<m; i++) 
	for (j=0; j<d; j++) {

	  tz.mul_2si(C(i,j),def);
	  Af(i,j).getData()=tz.get_d();  
	  
	}

      for (i=m; i<m+d; i++) 
	for (j=0; j<d; j++)
	  Af(i,j).getData()=(C(i,j)).get_d();  
	  
     {
	//ICI

	FP_NR<double> td;

	ZZ_mat<mpz_t> Afz;
	Afz.resize(d+1,d);
	
	for (i=0; i<m+d; i++) 
	  for (j=0; j<d; j++) {
	    td=(Af(i,j)).get_d();
	    (Afz(i,j)).set_f(td);
	  }

	double t,u,v,w;
	ratio(Afz,t,u,v,w);

	cout << endl << "---------------" << endl;
	
	cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	cout << ".. Average diagonal ratio: " << u << endl;
	cout << ".. Max diagonal ratio: " << v << endl;



      }

    
      B.assign(Af);

      B.hlll(delta);

      Af=B.getbase();
      
{
	//ICI

	FP_NR<double> td;

	ZZ_mat<mpz_t> Afz;
	Afz.resize(d+1,d);
	
	for (i=0; i<m+d; i++) 
	  for (j=0; j<d; j++) {
	    td=(Af(i,j)).get_d();
	    (Afz(i,j)).set_f(td);
	  }

	double t,u,v,w;
	ratio(Afz,t,u,v,w);

	cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	cout << ".. Average diagonal ratio: " << u << endl;
	cout << ".. Max diagonal ratio: " << v << endl;



      }
      
      Uf = B.getU();

      for (i=0; i<d; i++) 
	for (j=0; j<d; j++) {
	  tf = Uf(i,j).getData(); 
	  U(i,j).set_f(tf);  // Pour long double ou autre, vérifier et passer par set_z ? 
	}


      matprod_in_si(C,U);

      int size_of_U = maxbitsize(U,0,d,d);
      cout << "size of U: " << size_of_U << endl; 

      
      //print2maple(C,d+1,d);

      
    } // End Lehmer loop on the defect 
      
  } // End else lsigma > 0 
    
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
relations_lll_mixed(ZZ_mat<ZT>& C, matrix<FP_NR<mpfr_t> > F, int prec, int lalpha, int lsigma) { 

  
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

   epsilon=B.shift_epsilon(10);
   //epsilon=B.shift_epsilon(1000);  // Pour B4 

   cout << endl << "shifted epsilon = " << epsilon << endl << endl;

   FP_NR<mpfr_t> tmpabs; 

   // Unique shift 
   // ------------

   // ***** Régler le epsilon par rapport à la taille de U ?
   if (lsigma ==0) {  

     print2maple(B,n+1,n);
     
     Lattice<ZT, FT, MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> >, MatrixFT> D(B,NO_TRANSFORM,DEF_REDUCTION);
       
     D.hlll(0.99);
       
       B=D.getmixedbase();
	
       B.shift(-lalpha);

        print2maple(B,n+1,n);
              
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
     // Arrêter selon la précision !!!!!!!!!! et epsilon taille U !!!!!!!!

     // fplll
     ZZ_mat<ZT> A;
     
     A.resize(d+n,n);
       
     while (notfound) {
     //for (int K=0; K<1; K++) { 

       BL.set(B);

       current_shift+=lsigma;
       // ICI       
       cout << "current shift " << current_shift << endl; 
       
       BL.shift(current_shift);

       cout << "********** " << endl;
       
       print2maple(BL,4,n);
       
       Bt.mixed_put(BL, lsigma+2*n); // Heuristic 

       print2maple(Bt.getmixedbase(),4,n);
       
       print2maple(A,4,n);
       
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


  




/* ***********************************************

   Generalized Lehmer LLL  

   m x n  input matrix 

   shift : #bits, elementary shift (step) size à la Lehmer 
   unique shift if lsigma=0;

    WITH TRUNCATION 

   ********************************************** */

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
Lehmer_lll(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int shift=0, double delta=0.99) { 

  // TO DO FPLLL 

  if (shift == 0) { 
    
   
    Lattice<ZT, FT, MatrixZT, MatrixFT> B(A,NO_TRANSFORM,DEF_REDUCTION); 

    B.hlll(delta);

    C=B.getbase();

  } 
  else { 

    int m,n;
    int i,j,k;

    m=A.getRows();
    n=A.getCols();

    //print2maple(A,m,n);

    // Shifts to perform for every row 
    vector<vector<int> > shifts;   
    shifts.resize(m);

    int rowbits;
    vector<int> nbshifts; 
    nbshifts.resize(m);
    for (i=0; i<m; i++) nbshifts[i]=0;

    int nb;

    int maxnbshifts=0;

    // Loop for the shifts computation 
    for (i=0; i<m; i++) {

      rowbits=size_in_bits(A(i,0));
      for (j=1; j<n; j++) {
	rowbits=max(rowbits,size_in_bits(A(i,j)))-1;  // -1 for invertibility at the beginning 
      }

      //cout << endl;

      //cout << "row " << i << ": " << rowbits << endl; 

      nbshifts[i]=rowbits/shift;

      if((rowbits%shift) > 0) {
	nbshifts[i]+=1;
	nb=nbshifts[i];

	shifts[i].resize(nb+1);

	shifts[i][0]=-rowbits;
	for (k=1; k<nb; k++)
	  shifts[i][k]=shift;
	shifts[i][nb]=rowbits%shift;
	
      }
      else {
	nb=nbshifts[i];
	shifts[i].resize(nb+1);
	
	shifts[i][0]=-rowbits;
	for (k=1; k<=nb; k++)
	  shifts[i][k]=shift;
      }

      maxnbshifts=max(maxnbshifts,nbshifts[i]);

      //cout << "nbshifts " << nbshifts[i] << endl;
      //cout << "shifts " << shifts[i] << endl;

    } // end shifts computation   

    // ----------------
    // Main Lehmer loop 
    // ----------------

    Lattice<ZT, FT, MatrixZT, MatrixFT> Ct(A,TRANSFORM,NO_LONG);// DEF_REDUCTION);
    
    int s; // max global nb current shift 

    C.resize(m,n);
    
    for (i=0; i<m; i++)
      for(j=0; j<n; j++) 
	C(i,j)=A(i,j);

    vector<int> current_shift;
    current_shift.resize(m);
    for (i=0; i<m; i++) current_shift[i]=shifts[i][0];

    ZZ_mat<ZT> U;
    U.resize(n,n);

    ZZ_mat<ZT> T;
    T.resize(m,n);

    ZZ_mat<ZT> TT;
    TT.resize(1,n);

    // For fplll 
    ZZ_mat<ZT> V;
    V.resize(n,n);

    ZZ_mat<ZT> AT;
    AT.resize(n,m);

    //cout << "********** " << nbshifts << endl; 
    
    for (int row = m-1; row >= 0; row--) {

     
      for (s=1; s <= nbshifts[row]; s++) {

	// ICI 
	int mmax=0;
	for (j=0; j<n; j++)
	   if (mmax <  size_in_bits(C(0,j))) mmax=size_in_bits(C(0,j));
	
	//	current_shift[row]+=shifts[row][s];
	current_shift[row]=shift-mmax;
	
	cout << endl <<  current_shift[row] << endl; 
	
	// WITHOUT truncation test 
	Ct.shift_assign(C, current_shift, shift);
	
	{
	  // ICI
	  double t,u,v,w;
	  ratio(Ct.getbase(),t,u,v,w);
	  
	  cout << endl << "---------------" << endl;
	
	  cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	  cout << ".. Average diagonal ratio: " << u << endl;
	  cout << ".. Max diagonal ratio: " << v << endl;
	  cout << endl;
	  

	}

	{
	  // ICI

	  ZZ_mat<ZT> D;
	  D.resize(m,n);

	  D=Ct.getbase();
	  
	  ZZ_mat<long> E;
	  E.resize(12,12);

	  for (int i=0; i<12; i++)
	    for (int j=0; j<12; j++)
	      E(i,j)=size_in_bits(D(i,j));

	   int mmax=0;
	  for (int i=0; i<n; i++)
	    if (mmax <  size_in_bits(D(0,i))) mmax=size_in_bits(D(0,i));

	  cout << "max: " << mmax; 
	  //cout << E << endl; 


	}

	// fplll
	// -----
	
	// transpose(AT,Ct.getbase());
	
	// for (i=0; i<n; i++) 
	//   for (j=0; j<n; j++) 
	//     V(i,j)=0;

	// for (i=0; i<n; i++) 
	//   V(i,i)=1;
	
	// lllReduction(AT, V, delta, 0.51, LM_FAST,FT_DEFAULT,0);
	
	// transpose(U,V);
	
	// matprod_in(C,U);
	
	// HPLLL
	// -----
	
	Ct.hlll(delta);
	  
	U = Ct.getU();

	matprod_in(C,U);

	{
	  // ICI

	  ZZ_mat<ZT> D;
	  D.resize(m,n);

	  D=Ct.getbase();
	  
	 

	  int mmax=0;
	  for (int i=0; i<n; i++)
	    if (mmax <  size_in_bits(D(0,i))) mmax=size_in_bits(D(0,i));

	  cout << "max red: " << mmax << endl; 
	 


	}

	{
	  // ICI

	  ZZ_mat<ZT> D;
	  D.resize(m,n);

	  D=C;
	  
	  ZZ_mat<long> E;
	  E.resize(12,12);

	  for (int i=0; i<12; i++)
	    for (int j=0; j<12; j++)
	      E(i,j)=size_in_bits(D(i,j));


	  int mmax=0;
	  for (int i=0; i<n; i++)
	    if (mmax <  size_in_bits(D(0,i))) mmax=size_in_bits(D(0,i));

	  cout << "max: " << mmax << endl;
	  cout << "delta: " << mmax+current_shift[row] << endl; 
	  //cout << endl << E << endl; 


	}


      } // end on s - main Lehmer loop on the shifts 
    
    } // end loop on rows 

    //C=Ct.getbase();   
  }

 
 
  return 0;

}



/* ***********************************************

   Lehmer like LLL  

   (d+m) x n  input matrix 

   Shifting the upper d rows part 

   lsigma : bits, elementary shift (step) size à la Lehmer 
   unique shift if lsigma=0;

    WITH TRUNCATION 

   ********************************************** */

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
lehmer_lll(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int d, int lsigma=0) { 


  if (lsigma==0) { 
  Lattice<ZT, FT, MatrixZT, MatrixFT> Bt(A,NO_TRANSFORM,DEF_REDUCTION); 

  Bt.hlll(0.99);

  C=Bt.getbase();

  } 
  else { 
    int n,m;
    int i=0,j=0;

    

    m=A.getRows()-d;
    n=A.getCols();
 
  
    int max2up; 
    max2up=maxbitsize(A,0,d,n);

    MatrixZT B;
    B.resize(d+m,n);

    for (i=0; i<d+m; i++)
      for(j=0; j<n; j++) 
	B.set(i,j,A(i,j));

    //print2maple(B,d+m,n);   

    int max2down;
    max2down=maxbitsize(A,d,m,n);

    int delta;
    delta=max2up-max2down; 


    vector<int> shifts;   
    int nbshifts=0;

    nbshifts=delta/lsigma;

    if((delta%lsigma) > 0) {
      nbshifts+=1;
      shifts.resize(nbshifts+1);

      shifts[0]=-delta;
      for (i=1; i<nbshifts; i++)
	shifts[i]=lsigma;
      shifts[nbshifts]=delta%lsigma;
     
    }
    else {
    
      shifts.resize(nbshifts+1);

      shifts[0]=-delta;
      for (i=1; i<=nbshifts; i++)
	shifts[i]=lsigma;
    }


    Lattice<ZT, FT, MatrixZT, MatrixFT> Bt(A,TRANSFORM,DEF_REDUCTION);
      
    int s;
    int current_shift=0;

    // Main Lehmer loop 
    // ----------------
    
    // ICI 
    //cout << "shifts " << shifts << endl; 

    ZZ_mat<ZT> BL;
    BL.resize(d+m,n);


    for (s=0; s<=nbshifts; s++)  {
      
      //cout << "-------------------- " << s << endl; 

      set(BL,B);

      // ICI 
      //cout << "------------- " << endl; 
      //cout << "Size of B " << maxbitsize(BL) << endl;

      current_shift+=shifts[s];
      cout << "current shift " << current_shift << endl; 

      //print2maple(BL,d+m,n);

      shift_in(BL, current_shift, d);

      //print2maple(BL,d+m,n);

      Bt.put(BL, d, lsigma+20); // Heuristic 

      //print2maple(Bt.getbase(),d+m,n);
      //print2maple(B,d+m,n);

      
      Bt.hlll(0.99);

      {
	// ICI
	double t,u,v,w;
	ratio(Bt.getbase(),t,u,v,w);

	
	cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	cout << ".. Average diagonal ratio: " << u << endl;
	cout << ".. Max diagonal ratio: " << v << endl;


      }
      
      // Required multiplication update when BL has been truncated 
      matprod_in(B,Bt.getU());
      
      // ICI 
      //cout << "Size of U " << maxbitsize(Bt.getU()) << endl;


      //ICI 
      //cout << "-------- " << endl; 
      //print2maple(Bt.getbase(), d+m,n);
      //print2maple(B, d+m,n);
      // if not truncated on should not multiply evrything 

    } // End Lehmer loop 

   
    C.resize(d+m,n);

    for (i=0; i<d+m; i++)
      for(j=0; j<n; j++) 
	C(i,j)=B.get(i,j);
      
  } // End else lsigma > 0 
    
  return 0;
}

// No perturbation if lcond < 0 
template<class ZT> int  
additive_perturbation(ZZ_mat<ZT>& B, ZZ_mat<ZT> C, int lcond) { 


  int m,n;
  int i,j,lnorm;

  Z_NR<ZT> norm;
  Z_NR<ZT> r;
  int rsize;

  m=C.getRows();
  n=C.getCols();
 
  B.resize(m,n); // Output 
  
  if (lcond < 0) {

    for(j=0; j<n; j++) 
      for (i=0; i<m; i++) {
	B(i,j)=C(i,j);

      }
  }
  else {

    for(j=0; j<n; j++)  {
    
      norm.mul(C(0,j),C(0,j));
      for (i=1; i<m; i++)
	norm.addmul(C(i,j),C(i,j));
	
      lnorm = size_in_bits(norm);
	
      rsize = lnorm/2-lcond+1; 

      if (rsize <= 0) {
	for (i=0; i<m; i++) 
	  B(i,j)=C(i,j);
	}
      else {
	for (i=0; i<m; i++) {
	  r.randb(rsize);
	  B(i,j).add(r,C(i,j));
	}
      }
    }

  } 

  return 0;
}

template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
newlehmer_lll(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int shift, double delta) { 


  int m,n;
  int i,j;

  m=A.getRows();
  n=A.getCols();
 
  C.resize(m,n);    // Basis that is being reduced 
  
  for (i=0; i<m; i++)
    for(j=0; j<n; j++) 
      C(i,j)=A(i,j);
  

  Lattice<ZT, FT, MatrixZT, MatrixFT> Bt(A,TRANSFORM,DEF_REDUCTION);

  ZZ_mat<ZT> B;   // Temporary basis that is perturb for local reduction 
  B.resize(m,n);

  Z_NR<ZT> r;


  // Descending Lehmer loop 
  // ----------------------

  int bb;
  int nblov_tot=0;

  bb=maxbitsize(C);

  for (int K=0; K < bb/shift+2; K++) {
  

    // ICI 
    cout << "---------- " << endl; 
 
    //additive_perturbation<ZT>(B,C,lcond);

    /*for (i=0; i<m; i++) 
      for (j=0; j<n; j++) {
	r.randb(bb-K*shift);
	B(i,j).add(r,C(i,j));
	}*/

    for (i=0; i<m; i++) 
      for (j=0; j<n; j++) 
	  B(i,j)=C(i,j);

    if ((bb-K*shift) > 0) 
      for (i=1; i<m; i++) {
	r.randb(bb-K*shift);
	B(i,i-1).add(r,B(i,i-1));
      }

    cout << "perturb bits: " <<  bb-K*shift << endl; 

    //print2maple(C,10,10);

    //print2maple(C,m,n);
    //print2maple(B,m,n);

    //ICI 
    //Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Btest1(B,NO_TRANSFORM,DEF_REDUCTION);
    //Btest1.cond();

    Bt.assign(B);
  
    Bt.hlll(delta);
    cout << "nblov " << Bt.nblov << endl;
    nblov_tot+= Bt.nblov;
    cout << "nblov_tot " << nblov_tot << endl;

    matprod_in(C,Bt.getU());

    int sc;
    sc=maxbitsize(C);
    cout << "size C " << sc  << endl;
    
    //print2maple(Bt.getU(),n,n);


    //bb=maxbitsize(C) ;

    //cout << "old bb -bb " << oldbb-bb << endl;
    //cout << "bit size C " << bb << endl;

  
   
    //ICI 
    //Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Btest2(C,NO_TRANSFORM,DEF_REDUCTION);
    //cc=Btest2.cond();

    //cout << oldcc << "   " << cc << endl; 
    
    //if (Bt.nblov == n-1) K=KK; // Break  

    //print2maple(C,m,n);

  } // end Lehmer loop 


  return 0;
}



#endif 
