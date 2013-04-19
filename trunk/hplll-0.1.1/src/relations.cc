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


#include "relations.h"

#ifndef HPLLL_RELATIONS_CC
#define HPLLL_RELATIONS_CC


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



  Fgas<RT, mpz_t, FT,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, MatrixFT > F(B,setprec,d);
  //Fgas<RT, mpz_t, dpe_t,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, matrixexp<double, dpe_t> > F(B,setprec,d);
  // F.shift_epsilon(20);
  // F.setprec_internal(100);

  //Fgas<RT, mpz_t, ldpe_t,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, matrixexp<long double, ldpe_t> > F(B,setprec);
  //Fgas<RT, mpz_t, mpfr_t,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > F(B,setprec);

  int flag_decomp;
  flag_decomp=F.decomp(1.4142, nbrel, LONG_MAX);

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



#endif 
