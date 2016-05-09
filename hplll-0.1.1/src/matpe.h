/* matrix plus exponent class 

Created Jeu 20 jan 2011 17:16:18 CET
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


#ifndef HPLLL_MATPE_H
#define HPLLL_MATPE_H

#include "defs.h"
#include "mixed-col.h" 

//#include "matmixed.h" 

namespace hplll {

// ******************************************************************
// Begin matrix exp class  
// ******************************************************************

template<class T> struct colexp{
    T* col; /* significand */
    int  exp; /* exponent */
  } ;


template<class T, class DPET>
class MatrixPE {

protected:
  int n, d;
  vector<vector<T> >M;


public:
  vector<int> exp;   // À remettre en privé (ici pour debug) 
  
  // An empty matrix 
 MatrixPE() : n(0), d(0) {}

  // rows x cols, all elements are initialized with the default constructor of T 
 MatrixPE(int rows, int cols) : n(0), d(0) {
    resize(rows, cols);
  }
  
  /** Sets the dimensions of this matrix, preserving as much as possible of the
      content. The value of new elements is undefined. */
  
  void clear() {
    n = d = 0;
    M.clear();
    exp.clear();
  }


  void resize(int rows, int cols) {

    int i,j;
    n = rows;
    d = cols;
    M.resize(d);

    exp.resize(d);

    for (j=0; j<d; j++) M[j].resize(n);      

    for (j=0; j<d; j++) 
      for (i=0; i<n; i++) 
	M[j][i]=0;

    for (j=0; j<d; j++) 
      exp[j]=0;

  }

  /** Returns the number of rows */
  int getRows() const {
    return n;
  }

  /** Returns the number of columns */
  int getCols() const {
    return d;
  }


  // Normalisation au sens dpe d'une colonne 

  inline void normalize(int j, int nmax);

  // Accès aux éléments 

  inline FP_NR<DPET>  get(int i, int j);

  inline FP_NR<DPET>  get_non_normalized(int i, int j); 

  // Juste un élement ponctuellement, attention (coût) change aligne toute la colonne 
  // à éviter donc sans doute ? 

  inline void  set(int i, int j, double x) { 

    int tmpexp,d; 
    int k;

    if (x==0)
      M[j][i]=0;   // Quel exposant ? 
    else {
      
      M[j][i]=frexp(x,&tmpexp);
      d=tmpexp-exp[j];

      if (d >=0) {
	for (k=0; k<i; k++) 
	  M[j][k]=ldexp(M[j][k],-d);
	for (k=i+1; k<n; k++) 
	  M[j][k]=ldexp(M[j][k],-d);
	exp[j]=tmpexp;
      }
      else 
	M[j][i]=ldexp(M[j][i],d);
    }
  };

#ifdef HPLLL_WITH_LONG_DOUBLE
  inline void  set(int i, int j, long double x) { 

    int tmpexp,d; 
    int k;

    if (x==0)
      M[j][i]=0;   // Quel exposant ? 
    else {
      
      M[j][i]=frexpl(x,&tmpexp);
      d=tmpexp-exp[j];

      if (d >=0) {
	for (k=0; k<i; k++) 
	  M[j][k]=ldexpl(M[j][k],-d);
	for (k=i+1; k<n; k++) 
	  M[j][k]=ldexpl(M[j][k],-d);
	exp[j]=tmpexp;
      }
      else 
	M[j][i]=ldexpl(M[j][i],d);
    }
  };
#endif 

  // Juste un élement ponctuellement, attention (coût) change aligne toute la colonne 
  // à éviter donc sans doute ? 

  inline void  set(int i, int j, FP_NR<DPET> xdpe); 


  // Mise à jour de la colonne j à partir d'un vecteur de Z_NR<mpz_t>
  // On suppose qu'on remplace toute la colonne avec le nouvel exposant 

  inline void setcol(int j, Z_NR<long>* b, int beg, int length);

  inline void setcol(int j, Z_NR<double>* b, int beg, int length);
  
  inline void setcol(int j, Z_NR<mpz_t>* b, int beg, int length);

  inline void setcol(int j, FP_NR<mpfr_t>* f, int beg, int length);


  // Mise à jour de la colonne j à partir d'un vecteur de FP_NR<dpe_t>
  // Pas de normalisation 
  // On suppose qu'on remplace toute la colonne avec le nouvel exposant 

  inline void setcol(int j, const FP_NR<DPET>* vd, int length);     // double et long double 

  // From mixed matrices 
  inline void setcol(int j, const mixed_col<FP_NR<mpfr_t>, Z_NR<mpz_t> > c, int beg, int length);   

  // version const ? 

  const inline colexp<T> getcol(const int j, const int i=0) { 

    colexp<T> c;

    c.col=&M[j][i];
    c.exp=exp[j];

    return c;
  };


  // Recopie une colonne et écrase dans la nouvelle 

  inline void setcol(const int j, const colexp<T> wcol, const int nmax) {

    if (nmax !=0) {

      exp[j] = wcol.exp;
      T* w = wcol.col;
      T* v=&M[j][0];

      for (int i=0; i< nmax; i++) 
	v[i]=w[i];
    }
  };


  inline void setcol(const int j, const colexp<T> wcol, const int beg, const int nmax) {

    if (nmax !=0) {

      exp[j] = wcol.exp;
      T* w = wcol.col;
      T* v=&M[j][0];

      for (int i=0; i< nmax; i++) 
	v[i]=w[i];
    }
  };

  // Vector operation :  div, this col := w/a  à partir de la ligne k pour une longueur nmax 
  // on écrase tout le reste ds la colonne destination


  inline void div(const int j, const int k, const colexp<T> wcol, const FP_NR<DPET> a, const int nmax);


  // Vector operation :  fma, this col j := w - z*a  à partir de la ligne k 
  // Mise à jour de l'exposant en tenant compte du reste de la colonne en haut 

  // Apparemment pas besoin d'aligner par rapport au haut de la colonne ? 

  inline void fmasub(const int j, const int k,  
		     const colexp<T> wcol, const colexp<T> zcol, const FP_NR<DPET> a, const int nmax); 


  // Swaps col j and col k 
  inline void colswap(const int j, const int k)  {     // j < k 
 
    M[j].swap(M[k]);
    int tmpexp;
    tmpexp=exp[k];
    exp[k]=exp[j];
    exp[j]=tmpexp;

  };

  // en place de la matrice elle-même addcol col j := colj + col k 
  // On ne tient pas compte de ce qu'il y a par ailleurs dans la colonne cible 

  inline void addcol(const int j, const int k, const int l); 


  // subcol col j := colj - col k 

  inline void subcol(const int j, const int k, const int l); 

  
  // subcol col j := colj - a* col k 

  inline void submulcol(const int j, const int k,  const  FP_NR<DPET> a, const int l); 

}; // End matrix class 

//++++++++++++++++++++++++++++
// Spécialisations : 

 // subcol col j := colj - col k 

template<>  inline void MatrixPE<double,dpe_t>::subcol(const int j, const int k, const int l)  {    

    if (l != 0) {
      double* v=&M[j][0];
      double* w=&M[k][0];

      int d,i;
      d=exp[j]-exp[k];

      if (d >=0) {
	for (i=0; i< l; i++) 
	  v[i]-= ldexp(w[i],-d);
      }
      else {
	for (i=0; i< l; i++) 
	  v[i] = ldexp(v[i],d) - w[i];
	exp[j]=exp[k];
      }
    }
  };
  
#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline void MatrixPE<long double,ldpe_t>::subcol(const int j, const int k, const int l)  {    

    if (l != 0) {
      long double* v=&M[j][0];
      long double* w=&M[k][0];

      int d,i;
      d=exp[j]-exp[k];

      if (d >=0) {
	for (i=0; i< l; i++) 
	  v[i]-= std::ldexp(w[i],-d);
      }
      else {
	for (i=0; i< l; i++) 
	  v[i] = std::ldexp(v[i],d) - w[i];
	exp[j]=exp[k];
      }
    }
  };
#endif 

  // en place de la matrice elle-même addcol col j := colj + col k 
  // On ne tient pas compte de ce qu'il y a par ailleurs dans la colonne cible 

template<>  inline void MatrixPE<double,dpe_t>::addcol(const int j, const int k, const int l)  {     // j < k 
 
    if (l != 0) {


      double* v=&M[j][0];
      double* w=&M[k][0];

      int d,i;
      d=exp[j]-exp[k];

      if (d==0) {

	for (i=0; i< l; i++) 
	  v[i]+=w[i];

      }
      else if (d >0) {
	for (i=0; i< l; i++) 
	  v[i]+= ldexp(w[i],-d);
      }
      else {
	
	for (i=0; i< l; i++) 
	  v[i] = ldexp(v[i],d) + w[i];
 
	exp[j]=exp[k];
      }
    }

  };

#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline void MatrixPE<long double,ldpe_t>::addcol(const int j, const int k, const int l)  {     // j < k 
 
    if (l != 0) {

      long double* v=&M[j][0];
      long double* w=&M[k][0];

      int d,i;
      d=exp[j]-exp[k];

      if (d==0) {
	for (i=0; i< l; i++) 
	  v[i]+=w[i];
      }
      else if (d >=0) {
	for (i=0; i< l; i++) 
	  v[i]+= std::ldexp(w[i],-d);
      }
      else {
	for (i=0; i< l; i++) 
	  v[i] = std::ldexp(v[i],d) + w[i];
	exp[j]=exp[k];
      }
    }

  };
#endif 


  // Vector operation :  fma, this col j := w - z*a  à partir de la ligne k 
  // Mise à jour de l'exposant en tenant compte du reste de la colonne en haut 

  // Apparemment pas besoin d'aligner par rapport au haut de la colonne ? 
  // oui, changé Ven  5 oct 2012 09:53:05 CEST par exemple 
  // pour Householder complet (ds Rkept hlll pas besoin apparemment mais quand 
  // complet, oui 

template<>  inline void MatrixPE<double,dpe_t>::fmasub(const int j, const int k,  
		     const colexp<double> wcol, const colexp<double> zcol, const FP_NR<dpe_t> a, const int nmax) {

    if (nmax !=0) {

      int wexp = wcol.exp;
      double* w = wcol.col;

      int zexp = zcol.exp;
      double* z = zcol.col;

      double* v = &M[j][k];

      int d;
    
      int i; 

      double manta;     // double et long double 
      manta = DPE_MANT(a.getData());

      if ( manta == 0.0) {    // Vrai test cf add dpe   wexp > expprod + DPE_BITSIZE ? 
	exp[j]=wexp;
	for (i=0; i<nmax; i++)  
	  v[i]=w[i]; 
      }
      else { // Le scalaire a est non nul 

	int expprod = DPE_EXP(a.getData())+zexp;
	
	d=wexp-expprod;
  
	if (d >=0) {
	  manta=ldexp(manta,-d);
	  for (i=0; i<nmax; i++)  
	    v[i]=w[i] - manta*z[i];
	  exp[j]=wexp;
	}
	else {
	  double pp;     // double vs long double 
	  pp=pow(2.0,d);

	  for (i=0; i<nmax; i++)  {
	    //v[i]=ldexp(w[i],d)-manta*z[i];
	    v[i]=pp*w[i]-manta*z[i];
	  }
          // Attention haut de la colonne aussi, cas général 
          double* v = &M[j][0];
	  for (i=0; i<k; i++)   
	    v[i]=pp*v[i];
	  
	  exp[j]=expprod; 
	}
    
      }  // a != 0 
    } // nmax != 0
  };

#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline void MatrixPE<long double,ldpe_t>::fmasub(const int j, const int k,  
		     const colexp<long double> wcol, const colexp<long double> zcol, const FP_NR<ldpe_t> a, const int nmax) {

    if (nmax !=0) {

      int wexp = wcol.exp;
      long double* w = wcol.col;

      int zexp = zcol.exp;
      long double* z = zcol.col;

      long double* v = &M[j][k];

      int d;
    
      int i; 
 
      long double manta;     // double et long double 
      manta = LDPE_MANT(a.getData());

      if ( manta == 0.0) {    // Vrai test cf add dpe   wexp > expprod + DPE_BITSIZE ? 
	exp[j]=wexp;
	for (i=0; i<nmax; i++)  
	  v[i]=w[i]; 
      }
      else { // Le scalaire a est non nul 

	int expprod = LDPE_EXP(a.getData())+zexp;
	
	d=wexp-expprod;
  
	if (d >=0) {
	  manta=ldexpl(manta,-d);
	  for (i=0; i<nmax; i++)  
	    v[i]=w[i] - manta*z[i];
	  exp[j]=wexp;
	}
	else {
	  long double pp;     // double vs long double 
	  pp=std::pow(2.0,d);
	  for (i=0; i<nmax; i++)  
	    //v[i]=ldexp(w[i],d)-manta*z[i];
	    v[i]=pp*w[i]-manta*z[i];

	  // Attention haut de la colonne aussi, cas général 
	  long double* v = &M[j][0];
	  for (i=0; i<k; i++)   
	  v[i]=pp*v[i];

	    exp[j]=expprod; 
	}
    
      }  // a != 0 
    } // nmax != 0
 
  };
#endif 



  // Juste un élement ponctuellement, attention (coût) change aligne toute la colonne 
  // à éviter donc sans doute ? 

template<>  inline void MatrixPE<double,dpe_t>::set(int i, int j, FP_NR<dpe_t> xdpe) { 

    int tmpexp,d; 
    int k;

    double x;
    x=DPE_MANT(xdpe.getData());

    if (x==0.0)
      M[j][i]=0.0;
    else {
      
      M[j][i]=x;

      tmpexp=DPE_EXP(xdpe.getData());
      d=tmpexp-exp[j];

      if (d >=0) {
	for (k=0; k<i; k++) 
	  M[j][k]=ldexp(M[j][k],-d);
	for (k=i+1; k<n; k++) 
	  M[j][k]=ldexp(M[j][k],-d);
	exp[j]=tmpexp;
      }
      else 
	M[j][i]=ldexp(M[j][i],d);
    }
  };

#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline void MatrixPE<long double,ldpe_t>::set(int i, int j, FP_NR<ldpe_t> xdpe) { 

    int tmpexp,d; 
    int k;

    long double x;
    x=LDPE_MANT(xdpe.getData());

    if (x==0.0)
      M[j][i]=0.0;
    else {
      
      M[j][i]=x;

      tmpexp=LDPE_EXP(xdpe.getData());
      d=tmpexp-exp[j];

      if (d >=0) {
	for (k=0; k<i; k++) 
	  M[j][k]=ldexpl(M[j][k],-d);
	for (k=i+1; k<n; k++) 
	  M[j][k]=ldexpl(M[j][k],-d);
	exp[j]=tmpexp;
      }
      else 
	M[j][i]=ldexpl(M[j][i],d);
    }
  };
#endif 
  // Accès aux éléments 

template<>  inline  FP_NR<dpe_t>  MatrixPE<double,dpe_t>::get(int i, int j) { 

    FP_NR<dpe_t> x;
    DPE_MANT(x.getData()) = M[j][i];
    DPE_EXP(x.getData()) = exp[j];
    dpe_normalize (x.getData());
    
    return x;

  }

#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline  FP_NR<ldpe_t>  MatrixPE<long double, ldpe_t>::get(int i, int j) { 

    FP_NR<ldpe_t> x;
    LDPE_MANT(x.getData()) = M[j][i];
    LDPE_EXP(x.getData()) = exp[j];
    ldpe_normalize (x.getData());
    
    return x;

  }
#endif 

template<>  inline FP_NR<dpe_t>   MatrixPE<double, dpe_t>::get_non_normalized(int i, int j) {

    FP_NR<dpe_t> x;
    DPE_MANT(x.getData()) = M[j][i];
    DPE_EXP(x.getData()) = exp[j];
    
    return x;

  }


 
#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline FP_NR<ldpe_t>   MatrixPE<long double, ldpe_t>::get_non_normalized(int i, int j) { // voir double ou long double 

    FP_NR<ldpe_t> x;
    LDPE_MANT(x.getData()) = M[j][i];
    LDPE_EXP(x.getData()) = exp[j];
    //dpe_normalize (x.getData());
    
    return x;

  }
#endif 

  // subcol col j := colj - a* col k 

template<>  inline void MatrixPE<double,dpe_t>::submulcol(const int j, const int k,  const  FP_NR<dpe_t> a, const int l)  {     // j < k 
 
    if (l != 0) {

      double* w = &M[j][0]; 
      double* z = &M[k][0];

      int d; 
      int i; 

      double manta;     // double et long double 
      manta=DPE_MANT(a.getData());

      int expprod = DPE_EXP(a.getData())+exp[k];
      d=exp[j]-expprod;
  
      if (d >=0) {
	manta=ldexp(manta,-d);
	for (i=0; i<l; i++)  
	  w[i] -= manta*z[i]; 
      }
      else {
	for (i=0; i<l; i++)  
	  w[i]=ldexp(w[i],d)-manta*z[i];
	exp[j]=expprod; 
      }
    }
  };

#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline void MatrixPE<long double,ldpe_t>::submulcol(const int j, const int k,  const  FP_NR<ldpe_t> a, const int l)  {     // j < k 
 
    if (l != 0) {

      long double* w = &M[j][0]; 
      long double* z = &M[k][0];

      int d; 
      int i; 

      long double manta;     // double et long double 
      manta=LDPE_MANT(a.getData());

      int expprod = LDPE_EXP(a.getData())+exp[k];
      d=exp[j]-expprod;
  
      if (d >=0) {
	manta=ldexpl(manta,-d);
	for (i=0; i<l; i++)  
	  w[i] -= manta*z[i]; 
      }
      else {
	for (i=0; i<l; i++)  
	  w[i]=ldexpl(w[i],d)-manta*z[i];
	exp[j]=expprod; 
      }
    }
  };
#endif 


// Normalisation au sens dpe 

template<>  inline void MatrixPE<double,dpe_t>::normalize(int j, int nmax)  {  

      vector<int> tmpexp(nmax);
      int oldexp;

      double* v = &M[j][0];

      int d;
      int  e;

    if (nmax !=0) {

      // Calcul pour le nouvel exposant 
      oldexp=exp[j];
      d=0;

      for (int i=0; i<nmax; i++) {
	
	if (v[i]==0.0)  tmpexp[i]=INT_MIN;
	else {
	  v[i]=DPE_FREXP (v[i], &e); 
	  tmpexp[i]=oldexp+e;
	}
	d=max(d,tmpexp[i]);
      }
      
      for (int i=0; i<nmax; i++) {
	v[i]=ldexp(v[i],tmpexp[i]-d);
      }

      exp[j]=d;
    
    }
  
};

#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline void MatrixPE<long double,ldpe_t>::normalize(int j, int nmax)  {  

      vector<int> tmpexp(nmax);
      int oldexp;

      long double* v = &M[j][0];

      int d;
      int  e;

    if (nmax !=0) {

      // Calcul pour le nouvel exposant 
      oldexp=exp[j];
      d=0;

      for (int i=0; i<nmax; i++) {
	
	if (v[i]==0.0)  tmpexp[i]=INT_MIN;
	else {
	  v[i]=LDPE_FREXP (v[i], &e); 
	  tmpexp[i]=oldexp+e;
	}
	d=max(d,tmpexp[i]);
      }
      
      for (int i=0; i<nmax; i++) {
	v[i]=ldexpl(v[i],tmpexp[i]-d);
      }

      exp[j]=d;
    
    }
  
};
#endif 

// Mise à jour de la colonne j à partir d'un vecteur de Z_NR<mpz_t>
// On suppose qu'on remplace toute la colonne 

template<>  inline void MatrixPE<double,dpe_t>::setcol(int j, Z_NR<mpz_t>* b, int beg, int l) {  

  if (l !=0) {
    long maxexp=INT_MIN;
    vector<long> tmpexp(beg+l);

    for (int i=beg; i<beg+l; i++) {

      M[j][i]=mpz_get_d_2exp (&tmpexp[i], b[i].getData());

      maxexp=max(maxexp,tmpexp[i]);
    }

    exp[j]=maxexp;

    for (int i=beg; i<beg+l; i++) {
      M[j][i]=ldexp(M[j][i],tmpexp[i]-maxexp);
    }
  }
};

#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline void MatrixPE<long double,ldpe_t>::setcol(int j, Z_NR<mpz_t>* b, int beg, int l) {  

  if (l !=0) {

    mpfr_t xzf;
    mpfr_init2(xzf,64);  // À voir 

    int maxexp=INT_MIN;
    vector<long> tmpexp(beg+l);

    for (int i=beg; i<beg+l; i++) {

      mpfr_set_z(xzf, b[i].getData(),GMP_RNDN);

      M[j][i]=mpfr_get_ld_2exp (&tmpexp[i], xzf,GMP_RNDN);
  
      if (tmpexp[i] > maxexp) maxexp=tmpexp[i];

    }

    exp[j]=maxexp;

    for (int i=beg; i<beg+l; i++) {
      M[j][i]=ldexpl(M[j][i],tmpexp[i]-maxexp);
    }
  }
};
#endif 

// Mise à jour de la colonne j à partir d'un vecteur de Z_NR<long> et Z_NR<double>
// On suppose qu'on remplace toute la colonne
// Faire pour long double

// TODO better here simply cut and paste of above 
 
template<>  inline void MatrixPE<double,dpe_t>::setcol(int j, Z_NR<long>* b, int beg, int l) {  

  if (l !=0) {
    long maxexp=INT_MIN;
    vector<long> tmpexp(beg+l);

    for (int i=beg; i<beg+l; i++) {

      M[j][i]=b[i].get_d_2exp (&tmpexp[i]);

      maxexp=max(maxexp,tmpexp[i]);
    }

    exp[j]=maxexp;

    for (int i=beg; i<beg+l; i++) {
      M[j][i]=ldexp(M[j][i],tmpexp[i]-maxexp);
    }
  }
};

 template<>  inline void MatrixPE<double,dpe_t>::setcol(int j, Z_NR<double>* b, int beg, int l) {  

  if (l !=0) {
    long maxexp=INT_MIN;
    vector<long> tmpexp(beg+l);

    for (int i=beg; i<beg+l; i++) {

      M[j][i]=b[i].get_d_2exp (&tmpexp[i]);

      maxexp=max(maxexp,tmpexp[i]);
    }

    exp[j]=maxexp;

    for (int i=beg; i<beg+l; i++) {
      M[j][i]=ldexp(M[j][i],tmpexp[i]-maxexp);
    }
  }
};
 
// Mise à jour de la colonne j à partir d'un vecteur de FP_NR<mpfr_t>
// On suppose qu'on remplace toute la colonne 

template<>  inline void MatrixPE<double,dpe_t>::setcol(int j, FP_NR<mpfr_t>* b, int beg, int l) {  

  if (l !=0) {
    int maxexp=INT_MIN;
    vector<long> tmpexp(beg+l);

    for (int i=beg; i<beg+l; i++) {

      M[j][i]=mpfr_get_d_2exp (&tmpexp[i], b[i].getData(),GMP_RNDN);
      if (tmpexp[i] > maxexp) maxexp=tmpexp[i];

    }

    exp[j]=maxexp;

    for (int i=beg; i<beg+l; i++) {
      M[j][i]=ldexp(M[j][i],tmpexp[i]-maxexp);
    }
  }
};

#ifdef HPLLL_WITH_LONG_DOUBLE
template<>  inline void MatrixPE<long double, ldpe_t>::setcol(int j, FP_NR<mpfr_t>* b, int beg, int l) {  

  if (l !=0) {
    int maxexp=INT_MIN;
    vector<long> tmpexp(beg+l);

    for (int i=beg; i<beg+l; i++) {

      M[j][i]=mpfr_get_ld_2exp (&tmpexp[i], b[i].getData(),GMP_RNDN);
      if (tmpexp[i] > maxexp) maxexp=tmpexp[i];

    }

    exp[j]=maxexp;

    for (int i=beg; i<beg+l; i++) {
      M[j][i]=ldexp(M[j][i],tmpexp[i]-maxexp);
    }
  }
};
#endif 


  // Mise à jour de la colonne j à partir d'un vecteur de FP_NR<dpe_t>
  // Pas de normalisation 
  // On suppose qu'on remplace toute la colonne avec le nouvel exposant +++++++++++++

  
template<> inline void  MatrixPE<double,dpe_t>::setcol(int j, const FP_NR<dpe_t>* vd, int nmax) {

  double* v = &M[j][0];
  
  int d;
  
  if (nmax !=0) {

    
    double v_mant[nmax];
    DPE_EXP_T v_exp[nmax];
    

    // Calcul pour le nouvel exposant
    d=INT_MIN;
    
    for (int i=0; i<nmax; i++) {
      v_exp[i]=DPE_EXP(vd[i].getData());
      d=max(d,v_exp[i]);
      v_mant[i]=DPE_MANT(vd[i].getData());
    }
        
    for (int i=0; i<nmax; i++) {
      v[i]=ldexp(v_mant[i],v_exp[i]-d);
    }
    
    exp[j]=d;
    
  }  // Longueur !=0

};

// From mixed matrices 
template<>  inline void MatrixPE<double, dpe_t>::setcol(int j, const mixed_col<FP_NR<mpfr_t>, Z_NR<mpz_t> > c, int beg, int l) {

  int dR; 
  dR=c.dimRT;
  //int mZ;
  //mZ=c.dimZT;  // Could be used for some check on the value of l 

  FP_NR<mpfr_t>* vR;
  vR=c.colRT;

  Z_NR<mpz_t>* vZ;
  vZ=c.colZT;

  if (l !=0) {

    int maxexp=INT_MIN;
    vector<long> tmpexp(beg+l);

    // mpfr vector conversion 
    // ----------------------

    for (int i=beg; ((i<beg+l) && (i<dR)); i++) {

      M[j][i]=mpfr_get_d_2exp (&tmpexp[i], vR[i].getData(),GMP_RNDN);
      if (tmpexp[i] > maxexp) maxexp=tmpexp[i];

    }

    // mpz vector conversion 
    // ----------------------

    for (int i=0; i < l+beg-dR; i++) {

      M[j][i+dR]=mpz_get_d_2exp (&tmpexp[i+dR], vZ[i].getData());
      if (tmpexp[i+dR] > maxexp) maxexp=tmpexp[i+dR];

    }

    // Exponent alignment 
    // ------------------

    exp[j]=maxexp;
    for (int i=beg; i<beg+l; i++) {
      M[j][i]=ldexp(M[j][i],tmpexp[i]-maxexp);
    }

  } // nonzero length 

}; 

#ifdef HPLLL_WITH_LONG_DOUBLE
// From mixed matrices 
template<>  inline void MatrixPE<long double, ldpe_t>::setcol(int j, const mixed_col<FP_NR<mpfr_t>, Z_NR<mpz_t> > c, int beg, int l) {

  int dR; 
  dR=c.dimRT;
  //int mZ;
  //mZ=c.dimZT;

  FP_NR<mpfr_t>* vR;
  vR=c.colRT;

  Z_NR<mpz_t>* vZ;
  vZ=c.colZT;

  mpfr_t xzf;
  mpfr_init2(xzf,64);  // À voir 

  if (l !=0) {

    int maxexp=INT_MIN;
    vector<long> tmpexp(beg+l);

    // mpfr vector conversion 
    // ----------------------

    for (int i=beg; ((i<beg+l) && (i<dR)); i++) {

      M[j][i]=mpfr_get_ld_2exp (&tmpexp[i], vR[i].getData(),GMP_RNDN);
      if (tmpexp[i] > maxexp) maxexp=tmpexp[i];

    }


    // mpz vector conversion 
    // ----------------------

    for (int i=0; i < l+beg-dR; i++) {
      
      mpfr_set_z(xzf, vZ[i].getData(),GMP_RNDN);
      M[j][i+dR]=mpfr_get_ld_2exp (&tmpexp[i+dR], xzf,GMP_RNDN);

      if (tmpexp[i+dR] > maxexp) maxexp=tmpexp[i+dR];

    }


    // Exponent alignment 
    // ------------------

    exp[j]=maxexp;
    for (int i=beg; i<beg+l; i++) {
      M[j][i]=ldexp(M[j][i],tmpexp[i]-maxexp);
    }

  } // non zero length 

}; 
#endif 

#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void  MatrixPE<long double, ldpe_t>::setcol(int j, const FP_NR<ldpe_t>* vd, int nmax) {   

  long double* v = &M[j][0];

  int d;
  
  if (nmax !=0) {


    long double v_mant[nmax];
    LDPE_EXP_T v_exp[nmax];

    // Calcul pour le nouvel exposant
    d=INT_MIN;
    
    for (int i=0; i<nmax; i++) {
      v_exp[i]=LDPE_EXP(vd[i].getData());
      d=max(d,v_exp[i]);
      v_mant[i]=LDPE_MANT(vd[i].getData());
    }
        
    for (int i=0; i<nmax; i++) {
      v[i]=ldexpl(v_mant[i],v_exp[i]-d);
    }
    
    exp[j]=d;
    
    
  }  // Longueur !=0 

};
#endif 

  // Vector operation :  div, this col := w/a  à partir de la ligne k pour une loingueur nmax 
  // on écrase tout le reste ds la colonne destination

template<> inline void  MatrixPE<double,dpe_t>::div(const int j, const int k,  
				    const colexp<double> wcol, const FP_NR<dpe_t> a, const int nmax)  {

    if (nmax != 0) {
      int wexp = wcol.exp;
      double* w = wcol.col;

      double* v = &M[j][k];

      exp[j]=wexp-DPE_EXP(a.getData());

      const double b=1.0/DPE_MANT(a.getData());

      for (int i=0; i<nmax; i++)  
	v[i]=w[i]*b;

    }
  };

#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void  MatrixPE<long double, ldpe_t>::div(const int j, const int k,  
				    const colexp<long double> wcol, const FP_NR<ldpe_t> a, const int nmax)  {

    if (nmax != 0) {
      int wexp = wcol.exp;
      long double* w = wcol.col;

      long double* v = &M[j][k];

      exp[j]=wexp-LDPE_EXP(a.getData());

      const long double b=1.0/LDPE_MANT(a.getData());

      for (int i=0; i<nmax; i++)  
	v[i]=w[i]*b;

    }
  };

#endif 

// ******************************************************************
// matrix operations 
// ******************************************************************

// Norme_2, flottant cf srqt 

template<class T, class DPET> inline void fp_norm(FP_NR<DPET>& nn, const colexp<T> vcol, const int n);


template<> inline void fp_norm(FP_NR<dpe_t>& nn, const colexp<double> vcol, const int nmax)
{

  int vexp = vcol.exp;
  double* v = vcol.col;
  double tmpnn;

  if (nmax !=0) {

    DPE_EXP(nn.getData()) = 2*vexp;
    
    tmpnn=v[0]*v[0]; 
    for (int i=1; i<nmax; i++)  tmpnn+=v[i]*v[i];

    DPE_MANT(nn.getData())=tmpnn;
    dpe_normalize (nn.getData());

    dpe_sqrt(nn.getData(),nn.getData());

  }
}; 

#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void fp_norm(FP_NR<ldpe_t>& nn, const colexp<long double> vcol, const int nmax)
{

  int vexp = vcol.exp;
  long double* v = vcol.col;
  long double tmpnn;

  if (nmax !=0) {

    LDPE_EXP(nn.getData()) = 2*vexp;
    
    tmpnn=v[0]*v[0]; 
    for (int i=1; i<nmax; i++)  tmpnn+=v[i]*v[i];

    LDPE_MANT(nn.getData())=tmpnn;
    ldpe_normalize (nn.getData());

    ldpe_sqrt(nn.getData(),nn.getData());

  }
}; 
#endif 
 // Norme_2 au carré
 // Sans doute comme implanté ici pas très précis ? 

template<class T, class DPET> inline void fp_norm_sq(FP_NR<DPET>& nn, const colexp<T> vcol, const int n);


template<> inline void fp_norm_sq(FP_NR<dpe_t>& nn, const colexp<double> vcol, const int nmax) 
{

    int vexp = vcol.exp;
    double* v = vcol.col;
    double tmpnn;

    if (nmax !=0) {
      DPE_EXP(nn.getData()) = 2*vexp;
    
      tmpnn=v[0]*v[0]; 
      for (int i=1; i<nmax; i++)  {
           tmpnn+=v[i]*v[i];
          }

      DPE_MANT(nn.getData()) = tmpnn; 
      dpe_normalize (nn.getData());
    }
}; 

#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void fp_norm_sq(FP_NR<ldpe_t>& nn, const colexp<long double> vcol, const int nmax) 
{

    int vexp = vcol.exp;
    long double* v = vcol.col;
    long double tmpnn;

    if (nmax !=0) {
      LDPE_EXP(nn.getData()) = 2*vexp;
    
      tmpnn=v[0]*v[0]; 
      for (int i=1; i<nmax; i++)  tmpnn+=v[i]*v[i];

      LDPE_MANT(nn.getData()) = tmpnn; 
      ldpe_normalize (nn.getData());
    }
}; 
#endif 

// Scalar product 
template<class T, class DPET> inline void scalarprod(FP_NR<DPET>& nn, const colexp<T> vcol, const colexp<T> wcol, const int n);


template<> inline void scalarprod(FP_NR<dpe_t>& nn, const colexp<double> vcol, const colexp<double> wcol, const int nmax) { 

    double* v = vcol.col;
    double* w = wcol.col;

    double tmpnn;
    
    if (nmax !=0) {

      tmpnn=0.0;

      for (int i=0; i<nmax; i++) 
	tmpnn+=v[i]*w[i]; 

      DPE_EXP(nn.getData())=vcol.exp+wcol.exp;
      DPE_MANT(nn.getData()) = tmpnn; 
      dpe_normalize (nn.getData());
    }
};


#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void scalarprod(FP_NR<ldpe_t>& nn, const colexp<long double> vcol, const colexp<long double> wcol, const int nmax) { 

    long double* v = vcol.col;
    long double* w = wcol.col;

    long double tmpnn;

    if (nmax !=0) {

      tmpnn=0.0;

      for (int i=0; i<nmax; i++) 
	tmpnn+=v[i]*w[i]; 

      LDPE_EXP(nn.getData())=vcol.exp+wcol.exp;
      LDPE_MANT(nn.getData()) = tmpnn; 
      ldpe_normalize (nn.getData());
    }
};
#endif 
// Print in maple format 

template<class T, class DPET> void print2maple(MatrixPE<T,DPET> M, int n, int d) 
{
  
cout << "Matrix([";
 for (int i=0;i<n;i++) {
    cout << "[";
    for (int j=0;j<d-1;j++) {
      //(M.get(i,j)).print();
      hplllprint(M.get(i,j));
      cout << ", ";
    }
    //(M.get(i,d-1)).print();
    hplllprint(M.get(i,d-1));
    if (i<n-1) cout << "  ],\n"; else  cout << "  ]\n";
 }
  cout << "]);" << endl;

  };


template<class T, class DPET> void printcol(MatrixPE<T, DPET> M, int n, int beg, int k=0) 
{

cout << "Matrix([";
 for (int i=0;i<n;i++) {
    cout << "[";
    for (int j=beg;j<beg+k;j++) {
      M.get(i,j).print();
      cout << ", ";
    }
    M.get(i,beg+k).print();
    if (i<n-1) cout << "  ],\n"; else  cout << "  ]\n";
 }
  cout << "]);" << endl;

  };



 
void set_f(matrix<Z_NR<mpz_t> >& B, MatrixPE<double, dpe_t> R, long condbits)
{

  int n,d;

  n= B.getRows();
  d= B.getCols();

  FP_NR<dpe_t> norm,minval;
  minval=0.0;
  norm=0.0;

  int i,j;


  fp_norm(minval,R.getcol(0),n);

   // Avant Mar 29 avr 2014 10:42:12 CEST
  for (j=1; j<d; j++) {
    
    fp_norm(norm,R.getcol(j),n);
 
    if (minval.cmp(norm) > 0) minval=norm;
  }

 

  
  FP_NR<dpe_t> bf;

  Z_NR<mpz_t> tt,z;
  tt=0;
  z=0;
  FP_NR<dpe_t> rt;

  Z_NR<mpz_t> mm;
  mm=B(0,0);

  for (j=0; j<d; j++) {

    
    fp_norm(norm,R.getcol(j),n);
    
    for (i=0; i<n; i++) {

      bf.mul_2si(R.get(i,j),condbits+1);  // +1

      bf.div(bf,minval);
      B(i,j).set_f(bf);

      }
     
  }

}



#ifdef HPLLL_WITH_LONG_DOUBLE
void set_f(matrix<Z_NR<mpz_t> >& B, MatrixPE<long double, ldpe_t> R, long condbits)
{

  int n,d;

  n= B.getRows();
  d= B.getCols();

  FP_NR<ldpe_t> norm,minval;
  minval=0.0;
  norm=0.0;

  int i,j;


  fp_norm(minval,R.getcol(0),n);

   // Avant Mar 29 avr 2014 10:42:12 CEST
  for (j=1; j<d; j++) {
    
    fp_norm(norm,R.getcol(j),n);
 
    if (minval.cmp(norm) > 0) minval=norm;
  }


  FP_NR<ldpe_t> bf;

  Z_NR<mpz_t> tt,z;
  tt=0;
  z=0;
  FP_NR<ldpe_t> rt;

  Z_NR<mpz_t> mm;
  mm=B(0,0);

  for (j=0; j<d; j++) {

    
    fp_norm(norm,R.getcol(j),n);
    
    for (i=0; i<n; i++) {

      bf.mul_2si(R.get(i,j),condbits+1);  // +1
     
      bf.div(bf,minval);
      set_f(B(i,j),bf);

      }
     
  }

}
#endif


/*
template<class T> void set(MatrixPE<T>& B, matrix<T> A) 
{

  int m,n,i,j;

  m= A.getRows();
  n= A.getCols();

   for (i=0; i<m; i++) 
    for (j=0; j<n; j++) 
      B.set(i,j,A.get(i,j)); 
 
};
*/


} // end namespace hplll


#endif
