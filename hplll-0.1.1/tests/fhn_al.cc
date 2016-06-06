/* 
Ven  3 jui 2016 15:04:49 CEST
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


#include "hplll.h"

#include "wrappers.h"

#include <NTL/LLL.h>

using namespace NTL;

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {

  char results[]="benchmarks_results/fhn_al.results";    // ******** SPECIALIZE
  
  filebuf fb;
  iostream os(&fb);
  fb.open (results,ios::out);           


  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------

  int k,K;


  k=0;

  vector<int> d(100);

  k=0;

  //------------
 
  
  d[k]=128;
  k+=1;

 d[k]=200;
  k+=1;

 d[k]=256;
  k+=1;

  d[k]=300;
  k+=1;

  d[k]=340;
  k+=1;

  d[k]=380;
  k+=1;

  d[k]=400;
  k+=1;
 
  
  //-------------

  K=k;

  double delta=0.99;

  int run=0;
  
  Timer time;

  int status=0;

    os << endl << "(FPLLL,) HPLLL wrapper >= dim_prec_1, NTL running times / Alpha bases - long " << endl;   // ******** SPECIALIZE
    os << endl << "FPLL stopped at 200 (wrapper mpfr) " << endl;
    os << endl << "NTL FP infinite loop for 300 ==> XD" << endl; 
                                                                                    
    os <<         "-----------------------------------------------------------------------------" << endl << endl;
 
    for (int k=0; k<K; k++) { 

      int n=d[k];

      /*****************************************************************************/
      /*   i-th run  */
      /*****************************************************************************/
      
      run+=1;

      //--- Generation -----

      genalpha<mpz_t>(A,n,1.1);

      AT.resize(n,n);

      transpose(AT,A); 

      ZZ_mat<long> Along;
      matrix_cast(Along,A);

      // -------------------

      cout << n <<  endl; 
      
      cout << "--------------  HLLL" << endl << endl; 

      {

      	os << endl << "------------------------------------------------ " << endl ;


	Lattice<long, double, matrix<Z_NR<long> >,  matrix<FP_NR<double> > > B(Along,NO_TRANSFORM,DEF_REDUCTION);  //* name 

	time.start();

	verboseDepth=1;
	if (n <= DIM_PREC_1) status=B.hlll(delta); //* name
	
	else hlll<long>(tmpmat, A, 0.99, true, true);
      	verboseDepth=0;
	
	time.stop();

	 
      	os << "Run " << run << "  with dim = " << n << ",   delta = " << delta <<  endl << endl;
      	os << "    hlll: " << time << endl ;
      	time.print(os);
      	os << endl;

	if  (n <= DIM_PREC_1) matrix_cast(tmpmat,B.getbase());
	matrix_cast(tmpmat,B.getbase());
	
      	if (status ==0) {
      	  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(tmpmat,NO_TRANSFORM,DEF_REDUCTION); //* names

      	  T.isreduced(delta-0.1); //* name

	  double t,u,v,w;

	  ratio<mpz_t>(tmpmat,t,u,v,w);

	  cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	  cout << ".. Average diagonal ratio: " << u << endl;
	  cout << ".. Max diagonal ratio: " << v << endl;
	  cout << ".. First vector quality: " << w << endl;
      	} 
      	cout << endl; 

      	// cout << "--------------  FPLLL WRAPPER VERBOSE " << endl << endl; 

	if (n <=200) {
	  time.start();
	  lllReduction(AT, delta, 0.501, LM_WRAPPER,FT_DEFAULT,0,LLL_VERBOSE);
	  time.stop();
	  

	  os << "   fplll: " << time << endl << endl ;
	  time.print(os);
	  os << endl;
	  
	  transpose(tmpmat,AT);
	  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(tmpmat,NO_TRANSFORM,DEF_REDUCTION); //* name
	  T2.isreduced(delta-0.1); //* name
	}
	
	cout << "--------------  NTL  " << endl << endl; 

	Mat<ZZ> BN; 

	// Input basis 
	fb.close();
	fb.open ("tmp.txt",ios::out);
	os <<  transpose(A) ;
	fb.close();
	fb.open ("tmp.txt",ios::in);
	os >> BN;
	fb.close();
        system("rm tmp.txt");
	fb.open (results,ios::app);


	time.start();	

	if (n < 300) {
	  LLL_FP(BN,0.99,0,0,1); 
	}
	else
	  LLL_XD(BN,0.99,0,0,1);
	
	time.stop();

	fb.close();
	fb.open ("tmp.txt",ios::out);
	os <<  BN ;
	fb.close();
	fb.open ("tmp.txt",ios::in);
	os >> AT;
	fb.close();
	system("rm tmp.txt");
	fb.open (results,ios::app);



	transpose(tmpmat,AT);
  
	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(tmpmat,NO_TRANSFORM,DEF_REDUCTION);
	verboseDepth=0;
	T.isreduced(delta-0.1);


	os << "   ntl: " << time << endl << endl ;
      	time.print(os);
      	os << endl;
	
      } 
   
    }// End on runs, k loop


    // END 
    fb.close();


 
  return 0;
}
