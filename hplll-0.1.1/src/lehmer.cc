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
  lehmer_lll(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A, double delta, int lshift) { 
    
    verboseDepth-=1;
    
    if (lshift==0) {
      
      Lattice<mpz_t, FT, matrix<Z_NR<mpz_t> >, MatrixFT> B(A,NO_TRANSFORM,DEF_REDUCTION); 
      
    B.hlll(delta);
    
    C=B.getbase();
    
  }

  // Non trivial shift 
  // -----------------
  
  else {   // lsshift <> 0 

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

    //DBG
    print2maple(C,n,d);
    
    // Main Lehmer loop on the defect
    // ------------------------------
    
    while (def < 0) { 

      def = min(0,def+lshift);

     
      // Voir à optimiser la fonction
      // Faire pour 128 bits
      
      lift_truncate(Ct, C, def, 60); // Régler >= lshift 

      print2maple(Ct,n,d);

      B.assign(Ct);
      
      B.hlll(0.99);
 
      matprod_in_si(C,B.getU());

      print2maple(C,n,d);
      

      //DBG puis incrementer 
      //def=0;

      
//       // ICI 
//       int mmax=0;
//       for (j=0; j<d; j++)
// 	if (mmax <  size_in_bits(C(0,j))) mmax=size_in_bits(C(0,j));
	
//       def=lsigma-mmax;

//       cout << endl << "def: " << def << endl; 

	
      
//       for (i=0; i<m; i++) 
// 	for (j=0; j<d; j++) {

// 	  tz.mul_2si(C(i,j),def);
// 	  Af(i,j).getData()=tz.get_d();  
	  
// 	}

//       for (i=m; i<m+d; i++) 
// 	for (j=0; j<d; j++)
// 	  Af(i,j).getData()=(C(i,j)).get_d();  
	  
//      {
// 	//ICI

// 	FP_NR<double> td;

// 	ZZ_mat<mpz_t> Afz;
// 	Afz.resize(d+1,d);
	
// 	for (i=0; i<m+d; i++) 
// 	  for (j=0; j<d; j++) {
// 	    td=(Af(i,j)).get_d();
// 	    (Afz(i,j)).set_f(td);
// 	  }

// 	double t,u,v,w;
// 	ratio(Afz,t,u,v,w);

// 	cout << endl << "---------------" << endl;
	
// 	cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
// 	cout << ".. Average diagonal ratio: " << u << endl;
// 	cout << ".. Max diagonal ratio: " << v << endl;



//       }

    
//       B.assign(Af);

//       B.hlll(delta);

//       Af=B.getbase();
      
// {
// 	//ICI

// 	FP_NR<double> td;

// 	ZZ_mat<mpz_t> Afz;
// 	Afz.resize(d+1,d);
	
// 	for (i=0; i<m+d; i++) 
// 	  for (j=0; j<d; j++) {
// 	    td=(Af(i,j)).get_d();
// 	    (Afz(i,j)).set_f(td);
// 	  }

// 	double t,u,v,w;
// 	ratio(Afz,t,u,v,w);

// 	cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
// 	cout << ".. Average diagonal ratio: " << u << endl;
// 	cout << ".. Max diagonal ratio: " << v << endl;



//       }
      
//       Uf = B.getU();

//       for (i=0; i<d; i++) 
// 	for (j=0; j<d; j++) {
// 	  tf = Uf(i,j).getData(); 
// 	  U(i,j).set_f(tf);  // Pour long double ou autre, vérifier et passer par set_z ? 
// 	}


//       matprod_in_si(C,U);

//       int size_of_U = maxbitsize(U,0,d,d);
//       cout << "size of U: " << size_of_U << endl; 

      
//       //print2maple(C,d+1,d);

      
    } // End Lehmer loop on the defect 
      
  } // End else lshift > 0 
    
    verboseDepth+=1;
    return 0;
  }

} // end namespace hplll

#endif 
