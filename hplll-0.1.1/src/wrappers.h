/* Wrappers for reduction methods 

Created Lun 30 mai 2016 13:31:31 CEST  

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


#ifndef HPLLL_WRAPPERS_H
#define HPLLL_WRAPPERS_H

#include "hplll.h" 

namespace hplll { 


  // *******************************************************************************************
  
  // Wrapper for hlll, ZT used internally
  // ------------------------------------
  
  template<class ZT> Timer
    hlll(ZZ_mat<mpz_t>& C, const ZZ_mat<mpz_t> A, double delta = 0.99,
	 bool check=true, bool comput_cond=true);


  // *******************************************************************************************
  
  // Long int wrapper
  // ----------------
  
  template<> Timer
    hlll<long>(ZZ_mat<mpz_t>& C, const ZZ_mat<mpz_t> A, double delta,
	       bool check, bool comput_cond)
    {

      Timer time;

      time.start();
      
      verboseDepth-=1;
      
      typedef long ZT; // and cout below to change 
      
      int n=A.getRows();
      int d=A.getCols();
      

      int d1=DIM_PREC_1;
      int d2=DIM_PREC_2;
      
      // Reduction
      // *********
      
      ZZ_mat<ZT> B;
      B.resize(n,d);
      
      matrix_cast(B,A);  // From mpz_t to ZT 
      
      if (d <= d1) {

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, long int, double " << endl << endl; 
	}

	// Base LLL
	// --------
	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L(B,NO_TRANSFORM,DEF_REDUCTION);
      
	L.hlll(delta);
	
	matrix_cast(C,L.getbase());
      }
      else if (d <= d2) {

	ZZ_mat<ZT> T;

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, long int, double " << endl << endl; 
	}
      
	// Base LLL
	// --------
      
	T.resize(n,d1);

	set(T,B,n,d1);

	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L0(T,NO_TRANSFORM,DEF_REDUCTION);
       
	L0.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, long int, double -- Seysen reduction " << endl << endl; 
	}
	
	// Seysen
	// ------
      
	set(B,L0.getbase(),n,d1);

	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L(B,NO_TRANSFORM,SEYSEN_REDUCTION);
    

	L.hlll(delta);

	matrix_cast(C,L.getbase());

      }
      else {

	ZZ_mat<long> T;

      	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, long int, double " << endl << endl; 
	}
	
	// Base LLL
	// --------
      
	T.resize(n,d1);

	set(T,B,n,d1);

	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L0(T,NO_TRANSFORM,DEF_REDUCTION);
     
	L0.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, long int, double -- Seysen reduction " << endl << endl; 
	}
      
	// Seysen 
	// ------
      
	set(B,L0.getbase(),n,d1);

	T.resize(n,d2);

	set(T,B,n,d2);

	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L(T,NO_TRANSFORM,SEYSEN_REDUCTION);
     

	L.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, long int, long double -- Seysen reduction " << endl << endl; 
	}
      
	// Long double Seysen 
	// ------------------
      
	set(B,L.getbase(),n,d2);

	Lattice<ZT, long double, matrix<Z_NR<ZT> >,  matrix<FP_NR<long double> > > M(B,NO_TRANSFORM,SEYSEN_REDUCTION);
      
    
	M.hlll(delta);

	matrix_cast(C,M.getbase());
     

      }

      time.stop();
      
      // Reduction check
      // ---------------

      if (check) {

	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(C,NO_TRANSFORM,DEF_REDUCTION);

	L.isreduced(delta-0.1);

      }

      if (comput_cond) {

	double t,u,v,w;
	
	ratio(C,t,u,v,w);

	cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	cout << ".. Average diagonal ratio: " << u << endl;
	cout << ".. Max diagonal ratio: " << v << endl;
	cout << ".. First vector quality: " << w << endl << endl;
      }
      
      verboseDepth+=1;

      cout << endl << "Reduction time: " << time << endl << endl;
       
      return time;

  }

  // *******************************************************************************************
  
  // 128 int wrapper
  // ---------------
  
  template<> Timer
    hlll<__int128_t>(ZZ_mat<mpz_t>& C, const ZZ_mat<mpz_t> A, double delta,
	       bool check, bool comput_cond)
    {

      Timer time;

      time.start();
      
      verboseDepth-=1;
      
      typedef __int128_t ZT; // and cout below to change 
      
      int n=A.getRows();
      int d=A.getCols();

      int d1=DIM_PREC_1;
      int d2=DIM_PREC_2;

      // Reduction
      // *********
      
      ZZ_mat<ZT> B;
      B.resize(n,d);
      
      matrix_cast(B,A);  // From mpz_t to ZT 
      
      if (d <= d1) {

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, __int128_t, double " << endl << endl; 
	}

	// Base LLL
	// --------
	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L(B,NO_TRANSFORM,DEF_REDUCTION);
      
	L.hlll(delta);
	
	matrix_cast(C,L.getbase());
      }
      else if (d <= d2) {

	ZZ_mat<ZT> T;

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, __int128_t, double " << endl << endl; 
	}
      
	// Base LLL
	// --------
      
	T.resize(n,d1);

	set(T,B,n,d1);

	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L0(T,NO_TRANSFORM,DEF_REDUCTION);
       
	L0.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, __int128_t, double -- Seysen reduction " << endl << endl; 
	}
	
	// Seysen
	// ------
      
	set(B,L0.getbase(),n,d1);

	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L(B,NO_TRANSFORM,SEYSEN_REDUCTION);
    

	L.hlll(delta);

	matrix_cast(C,L.getbase());

      }
      else {

	ZZ_mat<ZT> T;

      	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, __int128_t, double " << endl << endl; 
	}
	
	// Base LLL
	// --------
      
	T.resize(n,d1);

	set(T,B,n,d1);

	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L0(T,NO_TRANSFORM,DEF_REDUCTION);
     
	L0.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl;
	  cout << "HLLL Wrapper, __int128_t, double -- Seysen reduction " << endl << endl; 
	}
      
	// Seysen 
	// ------
      
	set(B,L0.getbase(),n,d1);

	T.resize(n,d2);

	set(T,B,n,d2);

	Lattice<ZT, double, matrix<Z_NR<ZT> >,  matrix<FP_NR<double> > > L(T,NO_TRANSFORM,SEYSEN_REDUCTION);
     

	L.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "---------------------------------------------------------" << endl; 
	          cout << "HLLL Wrapper, __int128_t, long double -- Seysen reduction " << endl << endl; 
	}
      
	// Long double Seysen 
	// ------------------
      
	set(B,L.getbase(),n,d2);

	Lattice<ZT, long double, matrix<Z_NR<ZT> >,  matrix<FP_NR<long double> > > M(B,NO_TRANSFORM,SEYSEN_REDUCTION);
      
    
	M.hlll(delta);

	matrix_cast(C,M.getbase());
     

      }

      time.stop();
      
      // Reduction check
      // ---------------

      if (check) {

	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(C,NO_TRANSFORM,DEF_REDUCTION);

	L.isreduced(delta-0.1);

      }

      if (comput_cond) {

	double t,u,v,w;
	
	ratio(C,t,u,v,w);

	cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	cout << ".. Average diagonal ratio: " << u << endl;
	cout << ".. Max diagonal ratio: " << v << endl;
	cout << ".. First vector quality: " << w << endl << endl;
      }
      
      verboseDepth+=1;

      cout << endl << "Reduction time: " << time << endl << endl;
       
      return time;

  }


  // *******************************************************************************************
  
  // mpz_t wrapper 
  // -------------
  
  template<> Timer
    hlll<mpz_t>(ZZ_mat<mpz_t>& C, const ZZ_mat<mpz_t> A, double delta,
	       bool check, bool comput_cond)
    {

      Timer time;

      time.start();
      
      verboseDepth-=1;
      
      typedef mpz_t ZT; // and cout below to change 
      
      int n=A.getRows();
      int d=A.getCols();
      

      int d1=DIM_PREC_1;
      int d2=DIM_PREC_2;
            
      // Reduction
      // *********
      
     
      if (d <= d1) {

	if (verboseDepth >=0) {
	  cout << endl << "----------------------------------" << endl; 
	  cout << "HLLL Wrapper, mpz_t, double " << endl << endl; 
	}

	// Base LLL
	// --------
	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L(A,NO_TRANSFORM,DEF_REDUCTION);
      
	L.hlll(delta);
	
	C=L.getbase();
      }
      else if (d <= d2) {

	ZZ_mat<ZT> T;

	if (verboseDepth >=0) {
	  cout << endl << "----------------------------------" << endl; 
	  cout << "HLLL Wrapper, mpz_t, double " << endl << endl; 
	}
      
	// Base LLL
	// --------
	
	T.resize(n,d1);

	set(T,A,n,d1);

	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L0(T,NO_TRANSFORM,DEF_REDUCTION);
       
	L0.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "----------------------------------" << endl; 
	  cout << "HLLL Wrapper,  mpz_t, double -- Seysen reduction " << endl << endl; 
	}
	
	// Seysen
	// ------

	T.resize(n,d); 
      
	set(T,L0.getbase(),n,d1);
	
	for (int i=0; i<n; i++)
	  for (int j=d1; j<d; j++)
	    T(i,j)=A(i,j);

	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L(T,NO_TRANSFORM,SEYSEN_REDUCTION);
    
	L.hlll(delta);

	C=L.getbase();

	
      }
      else {

	ZZ_mat<ZT> T;

      	if (verboseDepth >=0) {
	  cout << endl << "----------------------------------" << endl; 
	  cout << "HLLL Wrapper,  mpz_t, double " << endl << endl; 
	}
	
	// Base LLL
	// --------
      
	T.resize(n,d1);

	set(T,A,n,d1);

	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L0(T,NO_TRANSFORM,DEF_REDUCTION);
     
	L0.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "----------------------------------" << endl; 
	  cout << "HLLL Wrapper,  mpz_t, double -- Seysen reduction " << endl << endl; 
	}
      
	// Seysen 
	// ------

	T.resize(n,d2); 
      
	set(T,L0.getbase(),n,d1);

	for (int i=0; i<n; i++)
	  for (int j=d1; j<d2; j++)
	    T(i,j)=A(i,j);

	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L(T,NO_TRANSFORM,SEYSEN_REDUCTION);
     

	L.hlll(delta);

	if (verboseDepth >=0) {
	  cout << endl << "----------------------------------" << endl; 
	  cout << "HLLL Wrapper,  mpz_t, long double -- Seysen reduction " << endl << endl; 
	}
      
	// Long double Seysen 
	// ------------------

	T.resize(n,d); 
      
	set(T,L.getbase(),n,d2);

	for (int i=0; i<n; i++)
	  for (int j=d2; j<d; j++)
	    T(i,j)=A(i,j);

	Lattice<mpz_t, ldpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<long double, ldpe_t> > M(T,NO_TRANSFORM,SEYSEN_REDUCTION);
      
   
	M.hlll(delta);

	C=M.getbase();
     

      }

      time.stop();
      
      // Reduction check
      // ---------------

      if (check) {

	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(C,NO_TRANSFORM,DEF_REDUCTION);

	L.isreduced(delta-0.1);

      }

      if (comput_cond) {

	double t,u,v,w;
	
	ratio(C,t,u,v,w);

	cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	cout << ".. Average diagonal ratio: " << u << endl;
	cout << ".. Max diagonal ratio: " << v << endl;
	cout << ".. First vector quality: " << w << endl << endl;
      }

      cout << endl << "Reduction time: " << time << endl << endl; 
      
      verboseDepth+=1;
      return time;

  }

  
} // end namespace hplll



#endif 
