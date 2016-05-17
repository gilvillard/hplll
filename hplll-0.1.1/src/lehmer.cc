/* Integer matrix nullspace  

Created Created Sam 11 mai 2013 15:52:10 CEST   
Copyright (C) 2013      Gilles Villard 
Restarted Mer 11 mai 2016 13:34:45 CEST

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



#include "hplll.h" 

#ifndef HPLLL_LEHMER_CC
#define HPLLL_LEHMER_CC


namespace hplll {



  /* ***********************************************

     Lehmer like LLL  

     n x d matrix, shift of the fisrt row  

     Under the Gaussian heuristic 

     ********************************************** */

  // Mettre LLL method avec fplll
  
  // reduction method avec Seysen 

  // Mettre un défaut au shift
  
  template<class ZT, class FT, class MatrixFT> int  
  lehmer_lll(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A, double delta, int lshift, bool truncate) { 
    
    Timer time,ttrunc,tred,tmul,ttot;
    
    verboseDepth-=1;
    
    if (verboseDepth >=0) {
      cout << endl << "Lehmer LLL " << endl << endl;  
    }
    
    // Trivial shift, direct LLL reduction
    // -----------------------------------
    
    if (lshift==0) {
      
      Lattice<mpz_t, FT, matrix<Z_NR<mpz_t> >, MatrixFT> B(A,NO_TRANSFORM,DEF_REDUCTION); 
      
    B.hlll(delta);
    
    C=B.getbase();
    
  }

  // Non trivial shift 
  // -----------------
  
  else {   // lshift <> 0 

   
    int i,j;

    int n=A.getRows();
    int d=A.getCols();

    
    C.resize(n,d); // Update of the basis 

    for (i=0; i<n; i++)
      for (j=0; j<d; j++)
	C(i,j)=A(i,j);

    ZZ_mat<ZT> Ct;  // Truncated basis  
    Ct.resize(n,d);  

    // Truncated lattice
    Lattice<ZT, FT, matrix<Z_NR<ZT> >, MatrixFT> B(Ct,TRANSFORM,DEF_REDUCTION);

    
    int bitsize = maxbitsize(A);

    int def = -bitsize;

    
    // Main Lehmer loop on the defect
    // ------------------------------
    
    while (def < 0) { 

      

      if (verboseDepth >=0) {
	cout << endl <<  "Current default: " << def << endl;  
      }
      
      // Truncation 
      // ----------
      time.start();

      def = min(0,def+lshift);
     
      // Voir à optimiser la fonction
      // Faire pour 128 bits

     
      
      if (truncate) 
	lift_truncate(Ct, C, def, 40); // Régler >= lshift 
      else
	lift_truncate(Ct, C, def, 0); //


       // DBG
      // {

      
      // 	double t,u,v,w;
      // 	ratio<mpz_t>(Ct,t,u,v,w);

      // 	cout << "log 2 Frobenius norm cond: " << t << endl;

      // 	//Lattice<ZT, mpfr_t,  matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > Ttest(T,NO_TRANSFORM,DEF_REDUCTION);
      // 	//Ttest.isreduced(delta-0.2);

      // }
      
      // DBG
      // {

      // 	cout << "max bit size: " << maxbitsize(Ct) << endl; 

      // }
      time.stop();
      ttot+=time;
      ttrunc+=time;


      // Reduction 
      // ---------
      time.start();

      B.assign(Ct);
      
      B.hlll(0.99);

      time.stop();
      ttot+=time;
      tred+=time;

      // DBG
      // {
      // cout << endl << endl << "*******************************" << endl << endl;

      // double t,u,v,w;
      // ratio<mpz_t>(B.getbase(),t,u,v,w);
      
      // cout << "log 2 Frobenius norm cond: " << t << endl;
      
      // //Lattice<ZT, mpfr_t,  matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > Ttest(T,NO_TRANSFORM,DEF_REDUCTION);
      // //Ttest.isreduced(delta-0.2);
      
      // }
     
     
      // Product
      // -------
      // Open mp à faire
      time.start();

      matprod_in_int(C,B.getU());
      
      
      time.stop();
      ttot+=time;
      tmul+=time;

    } // End Lehmer loop on the defect 
    
  } // End else lshift > 0 
    

    cout << endl;
    cout << "Total: " << ttot << endl;
    cout << "   liftings/truncations: " << ttrunc << endl;
    cout << "   reductions: " << tred << endl;
    cout << "   products: " << tmul << endl;
    
    verboseDepth+=1;
    return 0;
  }
  
} // end namespace hplll

#endif 
