/* matrix class 

Created Jeu 20 jan 2011 17:16:18   
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


#ifndef HPLLL_MAT_H
#define HPLLL_MAT_H

#include "defs.h"
#include "mixed-col.h"

// ******************************************************************
// Begin matrix class  
// ******************************************************************

// Comme un tableau lin�aire de colonnes 

template<class T>
class matrix {

protected:
  int n, d;
  vector<vector<T> >M;


public:
  // An empty matrix 
  matrix() : n(0), d(0) {}

  // rows x cols, all elements are initialized with the default constructor of T 
  matrix(int rows, int cols) : n(0), d(0) {
    resize(rows, cols);
  }

  /** Sets the dimensions of this matrix, preserving as much as possible of the
      content. The value of new elements is undefined. */

  void clear() {
    n = d = 0;
    M.clear();
  }


  // Met � z�ro
  void resize(int rows, int cols) {

    int i,j;
    n = rows;
    d = cols;
    M.resize(d);

    for (j=0; j<d; j++) M[j].resize(n);      

    for (j=0; j<d; j++) 
      for (i=0; i<n; i++) 
	M[j][i]=0;

  }

  /** Returns the number of rows */
  int getRows() const {
    return n;
  }

  /** Returns the number of columns */
  int getCols() const {
    return d;
  }

  // T& operator()(size_t i, size_t j): returns entry i,j, _rep[i*m+j]
  //- const T& operator()(size_t i, size_t j) const: same thing on a constant matrix 
  
  // a reference to the element (i, j). */
  inline T& operator()(int i, int j) {
    return M[j][i];
  }

  // 
  inline  T& get(const int i, const int j)                {return M[j][i];}

  // Dummy function for being generic with MatrixExp
  inline  T& get_non_normalized(const int i, const int j)                {return M[j][i];}

  // 
  inline  void set(const int i, const int j, const T a)                {M[j][i]=a;}

  // 
  inline  void set(const int i, const int j, const double a)                {M[j][i]=a;}

  // Returns a constant reference to the element (i, j) on constant objects. */
  inline const T& operator()(int i, int j) const {
    return M[j][i]; 
  }


  // 
  inline  T* getcol(const int j, const int i=0)                {return &M[j][i];}

  const inline T* getcol(const int j, const int i=0)  const   {return &M[j][i];}

  // Assumed to replace the whole column 
  inline  void setcol(const int j, const T* w, const int l)   {

    T* v=getcol(j);

    for (int i=0; i< l; i++) v[i]=w[i];
  }

// Assumed to replace the whole column 
  inline  void setcol(const int j, const T* w, const int beg, const int l)   {

    T* v=getcol(j);

    for (int i=beg; i< beg+l; i++) v[i]=w[i];
  }

  // Assumed to replace the whole column 
  inline void setcol(int j, Z_NR<mpz_t>* b, int beg, int l) {

    for (int i=beg; i<beg+l; i++) 
      set_z(M[j][i],b[i]);

 }

  // Assumed to replace the whole column 
  inline void setcol(int j, FP_NR<mpfr_t>* b, int beg, int l) {
    
    for (int i=beg; i<beg+l; i++) 
      set_mpfr(M[j][i],b[i]);
    
  }

  // Assumed to replace the whole column 
  // From mixed matrices 
  inline void setcol(int j, const mixed_col<FP_NR<mpfr_t>, Z_NR<mpz_t> > c, int beg, int l) {
    
    int dR; 
    dR=c.dimRT;
    //int mZ;   // Could be used for some check on the value of l
    //mZ=c.dimZT;
    
    FP_NR<mpfr_t>* vR;
    vR=c.colRT;
    
    Z_NR<mpz_t>* vZ;
    vZ=c.colZT;
    
    if (l !=0) {
            
      // mpfr vector conversion 
      for (int i=beg; ((i<beg+l) && (i<dR)); i++)
	set_mpfr(M[j][i],vR[i]);
      
      
      // mpz vector conversion 
      for (int i=0; i < l+beg-dR; i++) 
	set_z(M[j][i+dR],vZ[i]);
      
  } // non zero length 

}; 





  // Assumed to replace the whole column 
  inline void setcol(int j, FP_NR<double>* b, int beg, int l) {
    
    for (int i=beg; i<beg+l; i++) 
      (M[j][i]).set(b[i]);
    
  }

  // Assumed to replace the whole column 
  inline void setcol(int j, FP_NR<long double>* b, int beg, int l) {
    
    for (int i=beg; i<beg+l; i++) 
      (M[j][i]).set(b[i]);
    
  }

  // Assumed to replace the whole column 
  inline void setcol(int j, FP_NR<dpe_t>* b, int beg, int l) {
    
    for (int i=beg; i<beg+l; i++) 
      (M[j][i]).set(b[i]);
    
  }

  // Swaps col j and col k 
  inline void colswap(const int j, const int k)  {     // j < k 
 
    M[j].swap(M[k]);

  }

  // addcol col j := colj + col k 

  inline void addcol(const int j, const int k, const int l)  {     // j < k 
 
    T* v=getcol(j);
    T* w=getcol(k);

    for (int i=0; i< l; i++) 
      v[i].add( v[i], w[i]);
  }

  // subcol col j := colj - col k 

  inline void subcol(const int j, const int k, const int l)  {     // j < k 
 
    T* v=getcol(j);
    T* w=getcol(k);

    for (int i=0; i< l; i++) v[i].sub( v[i], w[i]);

  }
  
  // subcol col j := colj - a* col k 
  inline void submulcol(const int j, const int k,  const T a, const int l)  {     // j < k 
 
    T* v=getcol(j);
    T* w=getcol(k);

    for (int i=0; i< l; i++) v[i].submul(a,w[i]);

  }
  
// addmulcol col j := colj + a* col k 
  inline void addmulcol(const int j, const int k,  const T a, const int l)  {     // j < k 
 
    T* v=getcol(j);
    T* w=getcol(k);

    for (int i=0; i< l; i++) v[i].addmul(a,w[i]);

  }

  //ICI 29 mars 2013
  //   v := w - z*a
  inline void fmasub(const int j, const int k,  const T* w, const T* z, const T a, const int l)  {     // j < k 
    
    T* v=getcol(j)+k; // Check 

   
    for (int i=0; i<l; i++)  {
      
      v[i]=w[i];
      v[i].submul(z[i], a);    
    }

  }

  // Vector operation :  div, v := w/a  
  inline void div(const int j, const int k, T* w,  const T a, const int l) {

    T* v=getcol(j)+k;

    for (int i=0; i<l; i++)  v[i].div(w[i], a);

    };


  inline  void rowswap(const int i, const int k)   {
    
    T tmp;
    for (int j=0; j<d; j++) {
      tmp=M[j][i];
      M[j][i]=M[j][k];
      M[j][k]=tmp;

    }
  }

 inline  void addrow(const int i, const int k)   {
    
    for (int j=0; j<d; j++) 
      M[j][i].add(M[j][i],M[j][k]);
  }

inline  void subrow(const int i, const int k)   {
    
    for (int j=0; j<d; j++) 
      M[j][i].sub(M[j][i],M[j][k]);
  }

 inline  void addmulrow(const int i, const int k, const T a)   {
    
    for (int j=0; j<d; j++) 
      M[j][i].addmul(M[j][k],a);
  }


 inline  void shift(const int lalpha, const int nrows) { 
  
   for (int j=0; j<d; j++) 
     for (int i=0; i<nrows; i++) 
       M[j][i].mul_2si(M[j][i],lalpha);
    
}

inline void set(Matrix<T>& A) 
{

  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) 
      M[j][i]=A(i,j); 
 
}



}; // End matrix class 



// ******************************************************************
// matrix operations 
// ******************************************************************

// Norme_2, flottant cf srqt 
template<class T> inline void fp_norm(T&nn, const T* v, const int n) 
{

  nn.mul(v[0],v[0]); 

 for (int i=1; i<n; i++) 
   nn.addmul(v[i],v[i]);
  nn.sqrt(nn);

};

// Norme_2 au carré , flottant cf srqt 
// avec calcul de norme partielle pour le test de Lovasz 
template<class T> inline void fp_norm_sq(T&nn, const T* v, const int n) 
{

  nn.mul(v[0],v[0]); 
  for (int i=1; i<n; i++)  nn.addmul(v[i],v[i]);

};

// Scalar product 
template<class T> inline void scalarprod(T&nn, const T* v, const T* w, const int n) 
{

  nn.mul(v[0],w[0]); 

  for (int i=1; i<n; i++) 
    nn.addmul(v[i],w[i]);
};

// Vector operation :  addmul, v := v + a w  
template<class T> inline void vector_addmul(T* v, const T* w, const T a, const int n) 
{
  for (int i=0; i<n; i++)  v[i].addmul(w[i], a);
};

template<class T> inline void vector_addmul_ui(T* v, const T* w, const unsigned long int a, const int n) 
{
  for (int i=0; i<n; i++)  v[i].addmul_ui(w[i], a);
};

// Vector operation :  submul, v := v - a w  
template<class T> inline void vector_submul(T* v, const T* w, const T a, const int n) 
{
  for (int i=0; i<n; i++)  v[i].submul(w[i], a);

};

template<class T> inline void vector_submul_ui(T* v, const T* w, const unsigned long int a, const int n) 
{
  for (int i=0; i<n; i++)  v[i].submul_ui(w[i], a);

};

// Vector operation :  div, v := w/a  
template<class T> inline void vector_div(T* v,  const T* w, const T a, const int n) 
{
  for (int i=0; i<n; i++)  v[i].div(w[i], a);
};

// Vector operation :  mul, v := w*a  
template<class T> inline void vector_mul(T* v,  const T* w, const T a, const int n) 
{
  for (int i=0; i<n; i++)  v[i].mul(w[i], a);
};

// Vector operation :  sub, v := w -z   
template<class T> inline void vector_sub(T* v,  const T* w, const T* z, const int n) 
{
  for (int i=0; i<n; i++)  v[i].sub(w[i], z[i]);
};

// Vector operation :  fmasub, v := w - z*a   
template<class T> inline void vector_fmasub(T* v,  const T* w, const T* z, const T a, const int n) 
{

  for (int i=0; i<n; i++)  {
    v[i]=w[i];
    v[i].submul(z[i], a);    
  }


};


// ******************************************************************
// ******************************************************************

template<class T> inline  void shift_in(Matrix<T>& B, const int lalpha, const int nrows) { 
  
   for (int j=0; j<B.GetNumCols(); j++) 
     for (int i=0; i<nrows; i++) 
       (B(i,j)).mul_2si(B(i,j),lalpha);
      
}

template<class T> void print2maple(matrix<T> B, int n, int d) 
{
  // Stockage en lignes 
cout << "Matrix([";
 for (int i=0;i<n;i++) {
    cout << "[";
    for (int j=0;j<d-1;j++) {
      (B.get(i,j)).print();
      cout << ", ";
    }
    (B.get(i,d-1)).print();
    if (i<n-1) cout << "  ],\n"; else  cout << "  ]\n";
 }
  cout << "]);" << endl;

};

// Matrix assignment 
// -----------------

template<class T> void set(Matrix<T>& B, matrix<T> A) 
{

  int m,n,i,j;

  m= A.getRows();
  n= A.getCols();

   for (i=0; i<m; i++) 
    for (j=0; j<n; j++) 
      B(i,j)=A.get(i,j); 
 
};




// ******************************************************************
// Pas en place 
// ******************************************************************

template<class T> void matprod(Matrix<T>& C,  Matrix<T> B, Matrix<T> U) 
{

  int n,d,dres,i,j,k;

  n= B.GetNumRows();
  d= B.GetNumCols();
  dres=U.GetNumCols();

  for (i=0; i<n; i++)  {
    
    for (j=0; j<dres; j++) {

      C(i,j).mul(B(i,0),U(0,j));

      for (k=1; k<d; k++) {
	C(i,j).addmul(B(i,k),U(k,j));
      }
    }
  }

};


template<class T> void matprod_in(Matrix<T>& C, Matrix<T> U) 
{

  int m,n,i,j,k;

  m= C.GetNumRows();
  n= C.GetNumCols();

  Matrix<T> tmat;
  tmat.resize(m,n);

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) {
      tmat(i,j).mul(C(i,0),U(0,j));
      for (k=1; k<n; k++) {
	tmat(i,j).addmul(C(i,k),U(k,j));
      }
    }

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++)
      C(i,j)=tmat(i,j);
};


inline void matprod_in_si(ZZ_mat<mpz_t>& C, ZZ_mat<long int> U) 
{

  int m,n,i,j,k;

  m= C.GetNumRows();
  n= C.GetNumCols();

  Matrix<Z_NR<mpz_t> > tmat;
  tmat.resize(m,n);

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) {
      tmat(i,j).mul_si(C(i,0),U(0,j).GetData());
      for (k=1; k<n; k++) {
	tmat(i,j).addmul_si(C(i,k),U(k,j).GetData());
      }
    }

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++)
      C(i,j)=tmat(i,j);
};





template<class T> void matprod(matrix<T>& C,  matrix<T> B, matrix<T> U) 
{

  int n,d,dres,i,j,k;

  n= B.getRows();
  d= B.getCols();
  dres=U.getCols();

  for (i=0; i<n; i++)  {
    
    for (j=0; j<dres; j++) {

      C(i,j).mul(B(i,0),U(0,j));


      for (k=1; k<d; k++) {
	C(i,j).addmul(B(i,k),U(k,j));
      }
    }
  }

  };



template<class T> void matprod_in(matrix<T>& C, matrix<T> U) 
{

  int m,n,i,j,k;

  m= C.getRows();
  n= C.getCols();

  Matrix<T> tmat;
  tmat.resize(m,n);

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) {
      tmat(i,j).mul(C(i,0),U(0,j));
      for (k=1; k<n; k++) {
	tmat(i,j).addmul(C(i,k),U(k,j));
      }
    }

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++)
      C(i,j)=tmat(i,j);
};


template<class T> void matprod_in(matrix<T>& C, Matrix<T> U) 
{

  int m,n,i,j,k;

  m= C.getRows();
  n= C.getCols();

  Matrix<T> tmat;
  tmat.resize(m,n);

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) {
      tmat(i,j).mul(C(i,0),U(0,j));
      for (k=1; k<n; k++) {
	tmat(i,j).addmul(C(i,k),U(k,j));
      }
    }

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++)
      C(i,j)=tmat(i,j);
};


template<class T> void matprod(Matrix<T>& B, Matrix<T> U) 
{

  int n,d,i,j,k;

  n= B.GetNumRows();
  d= B.GetNumCols();

  Matrix<T> C;
  C.resize(n,d);

  for (i=0; i<n; i++)  {
    
    for (j=0; j<d; j++) {

      C(i,j).mul(B(i,0),U(0,j));


      for (k=1; k<d; k++) {
	C(i,j).addmul(B(i,k),U(k,j));
      }
    }
    for (j=0; j<d; j++) 
      B.Set(i,j,C(i,j));
  }

};

template<class T> void transpose(Matrix<T>& B, Matrix<T> A) 
{

  int m,n,i,j;

  m= A.GetNumRows();
  n= A.GetNumCols();

  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      B.Set(j,i,A(i,j));


};

template<class T> void print2maple(Matrix<T> B, int n, int d) 
{
  // Stockage en lignes 
cout << "Matrix([";
 for (int i=0;i<n;i++) {
    cout << "[";
    for (int j=0;j<d-1;j++) {
      B[i][j].print();
      cout << ", ";
    }
    B[i][d-1].print();
    if (i<n-1) cout << "  ],\n"; else  cout << "  ]\n";
 }
  cout << "]);" << endl;

};

// From fplll llldiff 

template <class T>
int matcmp (Matrix<T> B1, Matrix<T> B2, int c, int r)
{
  int test=1, i, j, sg;
  T tmp1;
  T tmp2;

  for (i=0; i<r; i++){
    sg = 1;
    tmp1.abs(B1.Get(i,0));
    tmp2.abs(B2.Get(i,0));
    if (tmp1.cmp(tmp2)!=0){
      //cerr << r << ", 0     " << tmp1 << "  " << tmp2 << "\n";
      test = 0;
    }
    if (tmp1.cmp(B1.Get(i,0))!=0) sg *=-1;
    if (tmp1.cmp(B2.Get(i,0))!=0) sg *=-1;

    if (sg == 1){
      for (j=1; j<c; j++){
        if (B1.Get(i,j).cmp(B2.Get(i,j))!=0){
          //cerr << i << " " << j << "     " << B1.Get(i,j) << "  " << B2.Get(i,j) << "\n";
          test = 0;
        }
      }
    }
    else{
      for (j=1; j<c; j++){
        tmp1.mul_si(B1.Get(i,j),-1);
        if (tmp1.cmp(B2.Get(i,j))!=0){
          //cerr << i << " " << j << "     " << B1.Get(i,j) << "  " << B2.Get(i,j) << "\n";
          test = 0;
        }
      }
    }
  }

  return (test);
}




template<class T> void print2maple(vector<vector<T> > B, int n, int d) 
{
  // Stockage en lignes 
cout << "Matrix([";
 for (int i=0;i<n;i++) {
    cout << "[";
    for (int j=0;j<d-1;j++) {
      B[i][j].print();
      cout << ", ";
    }
    B[i][d-1].print();
    if (i<n-1) cout << "  ],\n"; else  cout << "  ]\n";
 }
  cout << "]);" << endl;

};

template<class T> void printcol(matrix<T> B, int n, int beg, int k=0) 
{

cout << "Matrix([";
 for (int i=0;i<n;i++) {
    cout << "[";
    for (int j=beg;j<beg+k;j++) {
      B(i,j).print();
      cout << ", ";
    }
    B(i,beg+k).print();
    if (i<n-1) cout << "  ],\n"; else  cout << "  ]\n";
 }
  cout << "]);" << endl;

};

template<class T> void printcol(vector<vector<T> > B, int n, int beg, int k=0) 
{

cout << "Matrix([";
 for (int i=0;i<n;i++) {
    cout << "[";
    for (int j=beg;j<beg+k;j++) {
      B[i][j].print();
      cout << ", ";
    }
    B[i][beg+k].print();
    if (i<n-1) cout << "  ],\n"; else  cout << "  ]\n";
 }
  cout << "]);" << endl;

};

// ********************************************************************
// 
//       TRUNCATION OF A BASIS for lift, lehmer, nullspace, recursive 
//           n = m + d 
//           Upper part specific for knapsack or nullspace   
//           Lowerpart for the identity or the transformation matrix 
//
//  !!! Version  matrix<Z_NR<T> >
//      ==> Also modify ZZ_mat<T>  
//          Exactly identical function body
// 
// ********************************************************************       
// Check the roundings 
//
// Division by 2^sigma = 10 ^k the first m=n-d rows  
// And truncation, keeping t bits, 2^t = 10 ^ k 
// If t=0, keeps everything after the division 

/*
template<class T> void trunc_sigma(matrix<Z_NR<T> >& B, ZZ_mat<T> A, long n, long d, long t, long sigma);

template<> void trunc_sigma(matrix<Z_NR<mpz_t> >& B, ZZ_mat<mpz_t> A, long n, long d, long t, long sigma)

{

  int m=n-d;

  // Division of the first m rows  
  // and on the fly computation of the log[2] length 
  // ------------------------------------------------

  // Computation via pow10 = 10^k = 2^sigma 
  

  mpz_t pow10,ten;
  mpz_init(pow10); 
  mpz_init(ten);
  mpz_set_ui(ten, 10); 

  long exp10;
  double tmp;

  tmp=((double) sigma)*0.301029995;
  exp10=(long) tmp;
  mpz_pow_ui(pow10, ten, exp10);

  int i,j;

  // Et calcul de la longueur max en base 2 

  long max2=0;

  for (i=0; i<m; i++)
    for (j=0; j<d; j++)  {

      mpz_tdiv_q(B(i,j).getData(), A(i,j).getData(), pow10);

      max2=max(max2,(long) mpz_sizeinbase(B(i,j).getData(),2)); 
    }

  for (i=m; i<n; i++)  
    for (j=0; j<d; j++)
      max2=max(max2,(long) mpz_sizeinbase(A(i,j).getData(),2)); 


  // On tronque en en gardant msb t  
  // ------------------------------

  // To tune with the roundings, here keeps a bit more than t
  // If t=0, keeps everything, TODO: make it clean (do nothing in this case)

  // *****  Negative t
  // An identity part, corresponding transformation from scratch 
  if (t< 0) {

    for (i=m; i<n; i++)  
      for (j=0; j<d; j++)
	B(i,j)=0; 

    for (i=m; i<n; i++)  
      B(i,i-m)=1; 
  }
  // *****  Zero t
  // One keeps all the bits 
  else if (t==0) {

    for (i=m; i<n; i++)  
      for (j=0; j<d; j++)
	B(i,j)=A(i,j); 

  }
  // **** Positive t 
  // One keep t bits 
  // Unchanged w.r.t Lehmer 
  else {
  
    tmp=((double) max2-t)*0.301029995;
  
    exp10= max((long) 0, (long) tmp);
    mpz_pow_ui(pow10, ten, exp10);

    for (i=0; i<m; i++)
      for (j=0; j<d; j++)
	mpz_tdiv_q(B(i,j).getData(), B(i,j).getData(), pow10);

    for (i=m; i<n; i++)  
      for (j=0; j<d; j++)
	mpz_tdiv_q(B(i,j).getData(), A(i,j).getData(), pow10);
  }

  // Check to comment 
  //max2=0;
  //for (i=0; i<n; i++)  
    //for (j=0; j<d; j++)
      //max2=max(max2,(long) mpz_sizeinbase(B(i,j).getData(),2)); 

  //cout << "***** max bit after shift and truncation " << max2 << endl << endl; 
  // Compare max2 and t 

};
*/

// ********************************************************************
// 
//       TRUNCATION OF A BASIS for lift, lehmer, nullspace, recursive 
//           Upper part d rows  n columns 
//           Total number of rows m
//             hence lower part m-d rows 
//       
//       Throw tau bits in the upper part, by division bu 2^tau = 10^k
//         nothing thrown for tau = 0 
//       Then global truncation keeping t bits, 2^t = 10^k
//         if t=0 keeps everything
//
//      ==> Also modify matrix<Z_NR<T> >
//          Exactly identical function body
// 
// ********************************************************************       
// 
// Todo : check the roundings 
//

template<class ZT, class MatrixZT> void trunc(MatrixZT& B, ZZ_mat<ZT> A, long d, long n, long m, long t, long tau) {

  // Division of the first m rows  
  // and on the fly computation of the log[2] length 
  // ------------------------------------------------

  // Computation via pow10 = 10^k = 2^tau

  mpz_t pow10,ten;
  mpz_init(pow10); 
  mpz_init(ten);
  mpz_set_ui(ten, 10); 

  long exp10;
  double tmp;
  
  int i,j;

  // Division of the upper part 
  // and computation of the global max bit length 
 
  long max2=0;

 
  if (tau >0) {

    tmp=((double) tau)*0.301029995;
    exp10=(long) tmp;
    mpz_pow_ui(pow10, ten, exp10);
    
    for (i=0; i<d; i++)
      for (j=0; j<n; j++)  {
	mpz_tdiv_q(B(i,j).getData(), A(i,j).getData(), pow10);
	max2=max(max2,(long) mpz_sizeinbase(B(i,j).getData(),2));
      }
  }
  else {
    for (i=0; i<d; i++)
      for (j=0; j<n; j++) {  
	
	B.set(i,j,A(i,j)); 
	max2=max(max2,(long) mpz_sizeinbase(B(i,j).getData(),2));
      }
  }

  if (t >0) {

    for (i=d; i<m; i++)  
      for (j=0; j<n; j++)
	max2=max(max2,(long) mpz_sizeinbase(A(i,j).getData(),2)); 
  }

  // Truncation that keeps t bits  
  // ----------------------------

  // To tune with the roundings, here keeps a bit more than t
  // If t=0, keeps everything

  // *****  Negative t
  // An identity part, corresponding transformation from scratch 
  if (t< 0) {

    for (i=d; i<m; i++)  
      for (j=0; j<n; j++)
	B(i,j)=0; 

    for (i=d; i<m; i++)  
      B(i,i-d)=1; 
  }
  // *****  Zero t, one keeps evrything  
  else if ((t==0) || (t >= max2)) {

    for (i=d; i<m; i++)  
      for (j=0; j<n; j++)
	B(i,j)=A(i,j); 

  }
  // **** Positive t, one keep t bits 
  else {
  
    tmp=((double) max2-t)*0.301029995;
    exp10= max((long) 0, (long) tmp);
    mpz_pow_ui(pow10, ten, exp10);

    for (i=0; i<d; i++)
      for (j=0; j<n; j++)
	mpz_tdiv_q(B(i,j).getData(), B(i,j).getData(), pow10);

    for (i=d; i<m; i++)  
      for (j=0; j<n; j++)
	mpz_tdiv_q(B(i,j).getData(), A(i,j).getData(), pow10);
  }
  
};


/********************************************************/
/* ******         MAXBITSIZE        ******************* */
/********************************************************/

/*
template<class T> int maxbitsize(const ZZ_mat<T>& B);
*/

inline  int maxbitsize(const ZZ_mat<mpz_t>& B) {

  int l=0;

  int n=B.getRows();
  int d=B.getCols();

  for (int i=0; i<n ; i++) 
    for (int j=0; j<d; j++)
      l=max(l, (int) mpz_sizeinbase(B(i,j).getData(),2)); 

  return l;

}

template<class T> 
int maxbitsize(const ZZ_mat<T>& B, int d0, int d, int n) {

  int l=0;

  for (int i=d0; i<d ; i++) 
    for (int j=0; j<n; j++)
      l=max(l, size_in_bits(B(i,j))); 

  return l;

}


inline int
utime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec * 1000000 + rus.ru_utime.tv_usec ;
}

inline int
utimesec ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec ;
}

// Discover a triangular structure, integer matrix 
template<class T> inline void matrix_structure(vector<int>& structure, matrix<T> B, int n, int d)
{ 
  int i,k;
  structure.resize(d);

  // Bottom zeros 
  for (k=0; k<d; k++) {
    for (i=n-1; (i>=0) && (B(i,k)==0); i--) { } 
    structure[k]=i; 
  }
  // Make it triangular for limiting changes during the computation 
  for (k=1; k<d; k++)   structure[k]=max(structure[k-1],structure[k]);


}
 


#endif