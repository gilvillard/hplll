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


#ifndef HPLLL_LEHMER_CC
#define HPLLL_LEHMER_CC

// METTRE DANS LE NAMESPACE 
using namespace hplll; 

/* ***********************************************

   Lift LLL algorithm 


   ********************************************** */

// TO DO method lll
// Avec truncation
// Cas particuliers entre alpha et sigma 

template<class ZT, class FT, class MatrixFT> int  
relation_lift(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int alpha=0, double delta=0.99) { 

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

   int def = -bitsize;

   int target_def = -bitsize + alpha;
   
   int new_def;

   int found;

   FP_NR<FT> rel_bound;
   
   // Main loop on the shifts
   // -----------------------
   while (def < target_def) {

     int start;

     for (i=0; i<m; i++) 
       for (j=0; j<d; j++) 
	 L(i,j)=A_in(i,j);
    
    //current_shift += shift;

    start=utime();
         
     for (i=0; i<d; i++)
       for (j=0; j<d ; j++) 
     	 T(m+i,j).getData()=A_in(m+i,j).get_d(); // cf pb for assigning a double to Z_NR<double> 

     // Lattice<double, dpe_t,  matrix<Z_NR<double> >, MatrixPE<double, dpe_t> > Bp(T,TRANSFORM,DEF_REDUCTION,1);
     Lattice<double, FT,  matrix<Z_NR<double> >,  MatrixFT> Bp(T,TRANSFORM,DEF_REDUCTION,1);

     Bp.assignL(L);
     
     found = Bp.detect_lift(delta,def,target_def,new_def,rel_bound);

     start=utime()-start;
     cout << "   time A: " << start/1000 << " ms" << endl << endl;
     cout << " Def: " << new_def << endl;
     start=utime();
     
    
     
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
     
     start=utime()-start;
     cout << "   time B: " << start/1000 << " ms" << endl << endl;
     start=utime();

     cout << "***** " << new_def << endl;
     cout << "***** " << maxbitsize(A_in)+new_def << endl;
     
     cout << endl << "  size of U: " << maxbitsize(U,0,d,d)  << endl;
     

   }

   // found = 0
   cout << "**** There might not be relations of norm less than " << rel_bound << endl; 
  
   return 0;
  

} 


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

    Lattice<ZT, FT, MatrixZT, MatrixFT> Ct(A,TRANSFORM,DEF_REDUCTION,0);
    
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

    cout << "********** " << nbshifts << endl; 
    
    for (int row = m-1; row >= 0; row--) {

     
      for (s=1; s <= nbshifts[row]; s++) {

	current_shift[row]+=shifts[row][s];

	// WITHOUT truncation test 
	Ct.shift_assign(C, current_shift, shift);

	//print2maple(Ct.getbase(),m,n);
	
	transpose(AT,Ct.getbase());

	for (i=0; i<n; i++) 
	  for (j=0; j<n; j++) 
	    V(i,j)=0;

	for (i=0; i<n; i++) 
	  V(i,i)=1;

	lllReduction(AT, V, delta, 0.51, LM_FAST,FT_DEFAULT,0);

	transpose(U,V);

	matprod_in(C,U);
	
	//cout << "********************" << endl << endl;

	transpose(T,AT);
	//	print2maple(T,m,n);

	//print2maple(C,m,n);

	
	/* ** HPLLL 
	Ct.hlll(delta);
	  
	//cout << "Size of transform " << maxbitsize(Ct.getU()) << endl; 
	U = Ct.getU();

	matprod_in(C,U);
	*/

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
      //cout << "current shift " << current_shift << endl; 

      //print2maple(BL,d+m,n);

      shift_in(BL, current_shift, d);

      //print2maple(BL,d+m,n);

      Bt.put(BL, d, lsigma+8); // Heuristic 

      //print2maple(Bt.getbase(),d+m,n);
      //print2maple(B,d+m,n);
    
      Bt.hlll(0.99);

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
