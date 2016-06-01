/* Wrappers for slll without the Gaussian heuristic  

Created Lun 30 mai 2016 18:06:17 CEST 
Copyright (C) 2016      Gilles Villard 

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


#include "slll.h"


#ifndef HPLLL_SLLLWRAP_H
#define HPLLL_SLLLWRAP_H


namespace hplll { 

/* ***********************************************

        Wrapper for using slll, with gap 

   ********************************************** */


// Reduced until real column gap_status 
  
template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
slll_wrap_gap(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int gap_position, int S, double delta, int reduction_method=0) {

  int i,j;

  int gap_status;

  int n=A.getRows();
  int d=A.getCols();

  // First part of the basis
  // -----------------------
  
  ZZ_mat<ZT> B;

  B.resize(n,gap_position);

  for (i=0; i<n; i++)
    for (j=0; j<gap_position; j++)
      B.Set(i,j,A(i,j));
  
  // Reduction of the first part
  // ---------------------------
  
  SLattice<ZT, FT,  MatrixZT, MatrixFT>  LB(B,S,NO_TRANSFORM,reduction_method);

  gap_status=LB.hlll(delta,S,S,1000000);

  // Recursively
  // -----------
  if (gap_status >=2) {

    slll_wrap_gap<ZT, FT,  MatrixZT, MatrixFT>(B,LB.getbase(),gap_status,delta,reduction_method);
    
  }
  // For second part directly
  // ------------------------
  else {

    B=LB.getbase();
  }

  // Reduction of the second part
  // ----------------------------
  
  C.resize(n,d);

  for (j=0; j<gap_position; j++)
    for (i=0; i<n; i++)
      C.Set(i,j,B(i,j));

  for (j=gap_position; j<d; j++)
    for (i=0; i<n; i++)
      C.Set(i,j,A(i,j));
  
  Lattice<ZT, FT,  MatrixZT, MatrixFT>  LC(C,NO_TRANSFORM,reduction_method);
  
  LC.hlll(0.99);
  
  C=LC.getbase();

  return 0;
}





/* ***********************************************

        Wrapper for using slll  

   ********************************************** */


  
template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
slll_wrap(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int dthreshold, int S, double delta, int reduction_method=0) { 

  verboseDepth-=1;
  
  OMPTimer time,ttot,tinit,tsize,th,thlll;
  time.clear();
  ttot.clear();
  tinit.clear();
  tsize.clear();
  th.clear();
  thlll.clear();
  
  int i,j,k;

  int n=A.getRows();
  int d=A.getCols();

  ZZ_mat<ZT> B;

  int gap_status;
  
  int d1=min(d,dthreshold);

  // Initial reduction
  // -----------------
  
  B.resize(n,d1);

  for (i=0; i<n; i++)
    for (j=0; j<d1; j++)
      B.Set(i,j,A(i,j));

 
  Lattice<ZT, FT,  MatrixZT, MatrixFT>  L(B,NO_TRANSFORM,reduction_method);

  time.start();

  L.hlll(delta);
  
  time.stop();
  tinit+=time;
  
  C=L.getbase();

  cout << endl << "Initial  time: " << tinit << endl << endl;

  ttot+=tinit;
  
  // Rest of the reductions, one dimension at a time
  // -----------------------------------------------
  
  for (k = d1+1; k<=d; k++) {    // Real dims (not indices) 
      
    time.start();
    
    B.resize(n,k);
    
    for (i=0; i<n; i++)
      for (j=0; j<k-1; j++)
	B.Set(i,j,C(i,j));
    
    for (i=0; i<n; i++)
      B.Set(i,k-1,A(i,k-1));
    
    if (verboseDepth >= 0) 
      cout << "Discovering+ vector " << k  << "/" << d << endl;
    

    // Size reduction of the last column
    // ---------------------------------
    
    Lattice<ZT, FT,  MatrixZT, MatrixFT>  LR(B,NO_TRANSFORM,reduction_method);

    LR.householder_r(0);
    LR.householder_v(0);
    
    for (i=1; i<k-1; i++) {
      LR.householder_r(i);
      LR.householder_v(i);
    }

    if (reduction_method < 1) 
      LR.hsizereduce(k-1);
    else 
      LR.seysenreduce(k-1);
    
    time.stop();
    ttot+=time;

        
    if (verboseDepth >= 0) 
      cout << "     Size reduction: " << time << endl;


    // DBG 
    if (k >= 480)
      {
	cout << "dbg avant slll  " << k << endl;   

	cout << transpose(LR.getbase()) << endl;
	
      }
    

      
    time.start();

    
    // Reduction with the last column and recursive w.r.t. gap_status 
    // --------------------------------------------------------------

    
    SLattice<ZT, FT,  MatrixZT, MatrixFT>  L(LR.getbase(),S,NO_TRANSFORM,reduction_method);

    // DBG 
    if (k >= 480)
      {
	cout << "dbg apres slll  " << k << endl;   

	cout << transpose(L.getbase()) << endl;
	
      }
    
    gap_status=L.hlll(delta,S,S,1000000);
    
    if (gap_status >=2) {

      slll_wrap_gap<ZT, FT,  MatrixZT, MatrixFT>(C,L.getbase(),gap_status,S,delta,reduction_method);

    }

    else
      C=L.getbase();

    time.stop();

    ttot+=time;

    
    if (verboseDepth >=0) {
      cout << "     Phase+: " << time << endl;
      cout << "     Nblov: " << L.nblov << endl; 
      cout << "     Total: " << ttot << endl << endl;
    }

   
  } // End loop on k: extra columns 
  

  cout << endl; 
  cout << "Initial reduction: " << tinit << endl;
  cout << "Segment reduction: " << ttot << endl;

  verboseDepth+=1;
  return 0;

  }

 
} // end namespace hplll


#endif 


  
