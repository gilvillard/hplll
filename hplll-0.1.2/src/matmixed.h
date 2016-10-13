/* matrix exponent class 

Created Sam  4 mai 2013 15:15:46 CEST 
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

// ***********  PENSER CAS d=0 *********************

#ifndef HPLLL_MATRZ_H
#define HPLLL_MATRZ_H


#include "mat.h"
#include "mixed-col.h"

namespace hplll { 

template<template <class T> class MatrixT, class RT, class ZT>
  class MatrixRZ {
  
 protected:
  int d, m, n;
  MatrixT<RT> F;
  MatrixT<ZT> A;  
  
 public:
  
  // An empty matrix 
 MatrixRZ() : d(0), m(0), n(0) {}
  
 MatrixRZ(int rows1, int rows2, int cols) : d(0), m(0), n(0) {
    resize(rows1, rows2, cols);
  }


   /** Sets the dimensions of this matrix, preserving as much as possible of the
       content. The value of new elements is undefined. */

  void clear() {
    d = m = n = 0;
    F.clear();
    A.clear(); 
  }
  
  void resize(int rows1, int rows2, int cols) {
    
    d=rows1;
    m=rows2;
    n = cols;
    
    if (d>0) F.resize(d,n);
    A.resize(m,n);
    
  }

  /** Returns the number of rows */
  int getRowsRT() const {
    return d;
  }

   /** Returns the number of rows */
  int getRowsZT() const {
    return m;
  }
  
  /** Returns the number of columns */
  int getCols() const {
    return n;
  }
  
  inline  MatrixT<RT> getmatrixRT()  {return F;}
  
  inline  MatrixT<ZT> getmatrixZT()  {return A;}
  
  inline  RT& getRT(const int i, const int j)  {return F.get(i,j);}
  
  inline  ZT& getZT(const int i, const int j)  {return A.get(i,j);}
  
  // get for the integer part only 
  inline  ZT& get(const int i, const int j)  {return A.get(i,j);}


  inline  void set(const int i, const int j, const RT a) {
    F.set(i,j,a); 
  }

  inline  void setRT(const int i, const int j, const RT a) {
    F.set(i,j,a); 
  }
  
  inline  void set(const int i, const int j, const ZT a) {
    A.set(i-d,j,a);
  }
  
  inline  void setZT(const int i, const int j, const ZT a) {
    A.set(i,j,a); 
  }


//const inline mixed_col<RT,ZT> getcol(const int j, const int i=0) { 
  const inline mixed_col<RT,ZT> getcol(const int j) {
    
    mixed_col<RT,ZT> c;
    
    c.dimRT=d;
    c.dimZT=m;
    c.colRT=F.getcol(j); 
    c.colZT=A.getcol(j);
    
    return c;
  }


 const inline void  addcol(const int j, const int k, const int l) {
   
   if (d>0) {
     
     if (l <= d) {
       F.addcol(j,k,l);
     } 
     else {
       F.addcol(j,k,d);
       A.addcol(j,k,l-d);
     }
   }
   else 
     A.addcol(j,k,l);

 }


 const inline void  subcol(const int j, const int k, const int l) {
   
   if (d>0) { 
     if (l <= d) {
       F.subcol(j,k,l);
     } 
     else {
       F.subcol(j,k,d);
       A.subcol(j,k,l-d);
     }
   }
   else 
     A.subcol(j,k,l);
 }

 const inline void  submulcol(const int j, const int k, const ZT xz, const int l) {
   
   RT x;
   set_z(x,xz);
   
   if (d>0) {
     if (l <= d) {
       F.submulcol(j,k,x,l);
     } 
     else {
       F.submulcol(j,k,x,d);
       A.submulcol(j,k,xz,l-d);
     }
   }
   else 
     A.submulcol(j,k,xz,l);
 }

// GV Ven 11 oct 2013 15:48:27 CEST
const inline void  addmulcol_si(const int j, const int k, const long a, const int l) {
   
   RT x;
   Z_NR<long> xz;
   xz=a;

   set_z(x,xz);
   
   if (d>0) {
     if (l <= d) {
       F.addmulcol(j,k,x,l);
     } 
     else {
       F.addmulcol(j,k,x,d);
       A.addmulcol_si(j,k,a,l-d);
     }
   }
   else 
     A.addmulcol_si(j,k,a,l);
 }


//Ven 11 oct 2013 15:48:27 CEST
const inline void  addmulcol_si_2exp(const int j, const int k, const long a,  const long expo, const int l) {
   
   RT x;
   RT tmp;
   Z_NR<long> xz;
   xz=a;

   set_z(x,xz);
   
   if (d>0) {
     if (l <= d) {
       // prod par a puis prod par mul_2si 
       for (int i=0; i<l; i++) {
	 tmp=F.get(i,k);
	 tmp.mul(tmp,x);
	 tmp.mul_2si(tmp,expo);
	 tmp.add(F.get(i,j),tmp);
	 F.set(i,j,tmp);
       }
     } 
     else {
       for (int i=0; i<d; i++) {
	 tmp=F.get(i,k);
	 tmp.mul(tmp,x);
	 tmp.mul_2si(tmp,expo);
	 tmp.add(F.get(i,j),tmp);
	 F.set(i,j,tmp);
       }
      
       A.addmulcol_si_2exp(j,k,a,expo,l-d);
     }
   }
   else 
     A.addmulcol_si_2exp(j,k,a,expo,l);
 }





 const inline void  colswap(const int j, const int k) {
   
   if (d>0) F.colswap(j,k);
   A.colswap(j,k);
   
 }
 

 const inline void  set(MatrixRZ<MatrixT, RT, ZT> N) {
   
   for (int i=0; i<d; i++) 
     for (int j=0; j<n; j++) 
       F.set(i,j,N.getRT(i,j));
       
   for (int i=0; i<m; i++) 
     for (int j=0; j<n; j++) 
       A.set(i,j,N.getZT(i,j));

 } 


 // Multiply the RT part by 2 raised to lalpha
 void shift(const int lalpha) { 

  RT tt;
  
  for (int i=0; i<d; i++) 
    for (int j=0; j<n; j++) {
      tt.mul_2si(F.get(i,j),lalpha);
      F.set(i,j,tt);
    }
 };
 

 RT get_epsilon() {

   int setprec=(F.get(0,0)).getprec();
   
   int shiftepsilon=8;  // heuritic to check 
   
   mpfr_set_default_prec(100); // Temporary 
   RT r1, r2, ee, epsilon;
   double prec;
  
   prec=(double) shiftepsilon-setprec; 
   r1=prec;
   r2=2.0;
   r2.log(r2);
   r1.mul(r1,r2);
   ee.exponential(r1); 
   
   set_mpfr(epsilon,ee);
   mpfr_set_default_prec(setprec);
   
   return(epsilon);
 } 

RT shift_epsilon(int bits) {

   int setprec=(F.get(0,0)).getprec();
   
   int shiftepsilon=bits;  
   
   mpfr_set_default_prec(100); // Temporary 
   RT r1, r2, ee, epsilon;
   double prec;
  
   prec=(double) shiftepsilon-setprec; 
   r1=prec;
   r2=2.0;
   r2.log(r2);
   r1.mul(r1,r2);
   ee.exponential(r1); 
   
   set_mpfr(epsilon,ee);
   mpfr_set_default_prec(setprec);
   
   return(epsilon);
 } 

 
inline  int maxbitsizeRT() {

  int l=0;
  mpz_t data;
  mpz_init(data);

  for (int i=0; i<d ; i++) 
    for (int j=0; j<n; j++) {
      mpfr_get_z(data, (F.get(i,j)).get_data(), GMP_RNDN);
      l=max(l, (int) mpz_sizeinbase(data,2)); 
    }
  return l;

}

inline  int maxbitsizeZT() {

  int l=0;

  for (int i=0; i<m ; i++) 
    for (int j=0; j<n; j++)
      l=max(l, (int) mpz_sizeinbase((A.get(i,j)).get_data(),2)); 

  return l;

}


// ***********  PENSER CAS d=0 *********************

}; // ******** End mixed matrix class 


// Specializations 
// ***************



// Other functions 
// ***************


template<template <class T> class MatrixT, class RT, class ZT> void matprod_in(MatrixRZ<MatrixT,RT,ZT>& B, Matrix<ZT> U) {

  int d,m,n;
  int i,j;

  d=B.getRowsRT(); 
  m=B.getRowsZT();
  n=B.getCols();

  MatrixT<RT> F,UR;
  F.resize(d,n);
  UR.resize(n,n);

  RT tt;
  for (i=0; i< n; i++)
    for (j=0; j< n; j++) {
      tt.set_z(U(i,j));
      UR.set(i,j,tt);
    }

  matprod(F,B.getmatrixRT(),UR);

  for (i=0; i< d; i++)
    for (j=0; j< n; j++) {
      B.setRT(i,j,F.get(i,j));
    }

  MatrixT<ZT> A, UZ;
  A.resize(m,n);
  UZ.resize(n,n);

  for (i=0; i< n; i++)
    for (j=0; j< n; j++) 
      UZ.set(i,j,U(i,j));
    

  matprod(A,B.getmatrixZT(),UZ);

  for (i=0; i< m; i++)
    for (j=0; j< n; j++) {
      B.setZT(i,j,A.get(i,j));
    }

}

 
template<template <class T> class MatrixT, class RT, class ZT> 
  void mixed_trunc(MatrixRZ<MatrixT, RT,ZT>& B, MatrixRZ<MatrixT, RT,ZT> A, long t, long tau) {

  // Division of the first d rows  
  // and on the fly computation of the log[2] length 
  // ------------------------------------------------

  int d,m,n;

  d=A.getRowsRT(); 
  m=A.getRowsZT();
  n=A.getCols();


  // Computation via pow10 = 10^k = 2^tau

  mpz_t pow10,ten;
  mpz_init(pow10); 
  mpz_init(ten);
  mpz_set_ui(ten, 10); 

  long exp10;
  double tmp;
  

  ZT zpow10;
  
  
  RT mpow10;
  
  
  
  int i,j;
  
  // Division of the upper part 
  // and computation of the global max bit length 
  
  long max2=0;
  
  RT tt;


  if (tau >0) {

    tmp=((double) tau)*0.301029995;
    exp10=(long) tmp;
    mpz_pow_ui(pow10, ten, exp10);
    mpz_pow_ui(zpow10.get_data(), ten, exp10);
    set_z(mpow10,zpow10);

   
    for (i=0; i<d; i++)
      for (j=0; j<n; j++)  {
	tt.div(A.getRT(i,j), mpow10); 
	B.setRT(i,j,tt);
      }  
  } 
  else {
    for (i=0; i<d; i++)
      for (j=0; j<n; j++)  
	
	B.setRT(i,j,A.getRT(i,j));
  }
  
  max2=B.maxbitsizeRT();
  
  if (t >0) 
    max2=max((long) max2, (long) A.maxbitsizeZT()); 


  // Truncation that keeps t bits  
  // ----------------------------
  
  // To tune with the roundings, here keeps a bit more than t
  // If t=0, keeps everything

  // *****  Negative t
  // An identity part, corresponding transformation from scratch 
  if (t< 0) {
    
    ZT zero;
    zero=0;
    ZT one;
    one=1;
    
    for (i=0; i<m; i++)  
      for (j=0; j<n; j++)
	B.setZT(i,j,zero); 
    
    for (i=0; i<m; i++)  
      B.setZT(i,i,one); 
  }

  // *****  Zero t, one keeps evrything  
  else if (t==0) {
    
    for (i=0; i<m; i++)  
      for (j=0; j<n; j++)
	B.setZT(i,j,A.getZT(i,j)); 
  }
  // **** Positive t, one keep t bits
   
  else {
  
    tmp=((double) max2-t)*0.301029995;
    exp10= max((long) 0, (long) tmp);
    mpz_pow_ui(pow10, ten, exp10);

    mpz_pow_ui(zpow10.get_data(), ten, exp10);
    set_z(mpow10,zpow10);

    for (i=0; i<d; i++)
      for (j=0; j<n; j++) {
	tt.div(B.getRT(i,j), mpow10); 
	B.setRT(i,j,tt);
      }

    ZT tz,tzz;
    for (i=0; i<m; i++)  
      for (j=0; j<n; j++) {
	tz=A.getZT(i,j);
	mpz_tdiv_q(tzz.get_data(), tz.get_data(), pow10);
	B.setZT(i,j,tzz); 
      }
  }
  
}


  
template<template <class T> class MatrixT, class RT, class ZT> void print2maple(MatrixRZ<MatrixT, RT,ZT> B, int k, int l) 
{

  int d=B.getRowsRT();
 

  // Only real part 
  // --------------
  if (k < d) { 

    cout << "Matrix([";

    for (int i=0;i<k;i++) {
      cout << "[";
      for (int j=0;j<l-1;j++) {
	(B.getRT(i,j)).print();
	cout << ", ";
      }
      (B.getRT(i,l-1)).print();
      if (i<k-1) cout << "  ],\n"; else  cout << "  ]\n";
    }
   
    cout << "]);" << endl;
  } 
  // Both parts 
  else { 

    cout << "Matrix([";
    
    for (int i=0;i<d;i++) {
      cout << "[";
      for (int j=0;j<l-1;j++) {
	(B.getRT(i,j)).print();
	cout << ", ";
      }
      (B.getRT(i,l-1)).print();
      cout << "  ],\n"; 
    }

    for (int i=0;i<k-d;i++) {
      cout << "[";
      for (int j=0;j<l-1;j++) {
	(B.getZT(i,j)).print();
	cout << ", ";
      }
      (B.getZT(i,l-1)).print();
      if (i<k-d-1) cout << "  ],\n"; else  cout << "  ]\n";
    }
    
    cout << "]);" << endl;
    
  }
};



// Discover a triangular structure, integer matrix 
template<template <class T> class MatrixT, class RT, class ZT> inline 
  void matrix_structure(vector<int>& structure, MatrixRZ<MatrixT, RT, ZT> B, int d, int m, int n)
{ 
  int i,k;
  structure.resize(n);

 
  for (k=0; k<n; k++) {

    for (i=m-1; (i>=0) && (B.getZT(i,k)==0); i--) { } 
    structure[k]=i+d; 

    if (i == -1) {
      for (i=d-1; (i>=0) && ((B.getRT(i,k)).cmp(0.0) ==0); i--) { } 
      structure[k]=i; 
    }

  }

  // Make it triangular for limiting changes during the computation 
  for (k=1; k<n; k++)   structure[k]=max(structure[k-1],structure[k]);

}
 

} // end namespace hplll

#endif
