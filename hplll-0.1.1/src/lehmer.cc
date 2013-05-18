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



#include "hlll.cc" 


#ifndef HPLLL_LEHMER_CC
#define HPLLL_LEHMER_CC



/* ***********************************************

   Lehmer like LLL nullspace 

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
    
    cout << "shifts " << shifts << endl; 

    ZZ_mat<ZT> BL;
    BL.resize(d+m,n);


    for (s=0; s<=nbshifts; s++)  {
      
      cout << "---- " << s << endl; 

      set(BL,B);

      current_shift+=shifts[s];
      //cout << "current shift " << current_shift << endl; 

      shift_in(BL, current_shift, d);

      //print2maple(BL,d+m,n);

      Bt.put(BL, d, lsigma+n); // Heuristic 
     

      //print2maple(Bt.getbase(),d+m,n);
      //print2maple(B,d+m,n);
    
      Bt.hlll(0.99);

      // Required multiplication update when BL has been truncated 
      matprod_in(B,Bt.getU());
      // if not truncated on should not multiply evrything 

    } // End Lehmer loop 

   
    C.resize(d+m,n);

    for (i=0; i<d+m; i++)
      for(j=0; j<n; j++) 
	C(i,j)=B.get(i,j);
      
  } // End else lsigma > 0 
    
  return 0;
}

#endif 
