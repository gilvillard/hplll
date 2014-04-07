/* Very preliminary implementation of recursive L1, knapsack case 

Created Created   Mar 13 mar 2012 10:37:59 CEST   
Copyright (C) 2012, 2013      Gilles Villard 

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

#ifndef HPLLL_L1_H
#define HPLLL_L1_H

namespace hplll { 


// The transform, U, is allocated outside  

template<class ZT, class FT, class MatrixZT, class MatrixFT> 
void rec(ZZ_mat<mpz_t>& U, const ZZ_mat<mpz_t> B, Lattice<ZT, FT, MatrixZT, MatrixFT >& L, 	 
	 long leafbits, int l, int level, int lllmethod=HLLL) {

  int n,d,i,j;
  
  n=B.getRows();
  d=B.getCols();


  if (l <= leafbits) { 

     if (lllmethod == HLLL) { 
       
       L.put(B,1,0); 
       
       L.hlll(0.99);
       U=L.getU();
       
     } 
     else if (lllmethod == FPLLL) { 
       
       ZZ_mat<mpz_t> AT;
       AT.resize(d,n);

       for (i=0; i<d; i++) 
	 for (j=0; j<n; j++) 
	   AT(i,j)=B(j,i);
       
      
       ZZ_mat<mpz_t> V;
       V.resize(d,d);
       for (i=0; i<d; i++) 
	 V(i,i)=1;

       //lllReduction(AT, 0.99, 0.51, LM_HEURISTIC,FT_DPE,0);
       lllReduction(AT, V, 0.99, 0.51, LM_WRAPPER,FT_DEFAULT,0);
      
       for (i=0; i<d; i++) 
	 for (j=0; j<d; j++) 
	   U(i,j)=V(j,i);

     } 
     //print2maple(U,d,d);

  } // Endif leaf 

   
   else {
     ZZ_mat<mpz_t> C1,C2;

     C1.resize(n,d);
     C2.resize(n,d);

     //cout << "Call A " << l/2 << endl; 
     
    
     trunc<mpz_t, ZZ_mat<mpz_t> >(C1, B, 1, d, n, l/2+d+2*(level+1), l/2); 
     //trunc<mpz_t, ZZ_mat<mpz_t> >(C1, B, 1, d, n, 0, l/2); 
     

     rec(U, C1, L, leafbits, l/2, level+1, lllmethod);

     // ICI 
     /*if (level < 2) {
       cout << "*** " << maxbitsize(U) << "   " << maxbitsize(B) << endl; 
       }*/
     matprod(C1,B,U);
    

     ZZ_mat<mpz_t> U2;
     U2.resize(d,d);

     // cout << "Call B " << l/2 << endl; 
     trunc<mpz_t, ZZ_mat<mpz_t> >(C2, C1, 1, d, n, l/2+d+2*(level+1), 0); 
     
     rec(U2, C2, L, leafbits, l/2, level+1, lllmethod);

     //rec(U2, C1, L, l/2, level+1, lllmethod);
     
     
     matprod(U,U2); 
   

   }
  

};


 void l1(ZZ_mat<mpz_t>& B, const ZZ_mat<mpz_t> A, long leafbits, int lllmethod=HLLL) { 


   //Lattice<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > > L(A,TRANSFORM,DEF_REDUCTION);
   Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double,dpe_t> > L(A,TRANSFORM,DEF_REDUCTION);

   int n,d;
   int l;

   n=A.getRows();
   d=A.getCols();

   B.resize(n,d);

   ZZ_mat<mpz_t> U;
   U.resize(d,d);

   l=maxbitsize(A);

   rec(U, A, L, leafbits, l, 0, lllmethod);

   matprod(B,A,U);

};



} // end namespace hplll


#endif 
