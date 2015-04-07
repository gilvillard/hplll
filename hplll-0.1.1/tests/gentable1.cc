/* Integer matrix nullspace test file  

Created Dim  7 avr 2013 16:54:03 CEST
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


#include "hlll.h"
#include "matgen.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  filebuf fb;
  iostream os(&fb);
  fb.open ("table1.txt",ios::out);

  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------

  int k,K;

  vector<int> d(100);
  vector<int> bits(100);

  k=0;

  //------------
  d[k]=20;
  bits[k]=400;
  k+=1;

  d[k]=60;
  bits[k]=3600;
  k+=1;

  d[k]=80;
  bits[k]=4600;
  k+=1;
  
  // d[k]=120;
  // bits[k]=7200;
  // k+=1;

  // d[k]=160;
  // bits[k]=12800;
  // k+=1;

  // d[k]=180;
  // bits[k]=16000;
  // k+=1;

  // d[k]=220;
  // bits[k]=16000;
  // k+=1;

  // d[k]=240;
  // bits[k]=20000;
  // k+=1;

  // d[k]=260;
  // bits[k]=24000;
  // k+=1;

  //-------------

  K=k;

  double delta=0.99;

  int run=0;
  
  Timer time;


    os << endl << "FPLLL wrapper and HPLLL running times / intrel bases" << endl; 
    os <<         "----------------------------------------------------" << endl << endl;
 
    for (int k=0; k<K; k++) { 


      /*****************************************************************************/
      /*   i-th run  */
      /*****************************************************************************/
      
      run+=1;

      A.resize(d[k]+1,d[k]);
      AT.resize(d[k],d[k]+1);

      AT.gen_intrel(bits[k]);
      transpose(A,AT);


      cout << "--------------  HLLL" << endl << endl; 

      {
	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B1(A,NO_TRANSFORM,DEF_REDUCTION);  //* name 

	time.start();
	B1.hlll(delta); //* name
	time.stop();

	os << "Run " << run << "  with d = " << d[k] << ",  bits = " << bits[k] << ",  delta = " << delta <<  endl << endl;
	os << "    hlll: " << time << endl ;
    
	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T11(B1.getbase(),NO_TRANSFORM,DEF_REDUCTION); //* names

	T11.isreduced(delta-0.1); //* name

	cout << endl; 

	cout << "--------------  FPLLL WRAPPER" << endl << endl; 
    
	time.start();
	lllReduction(AT, delta, 0.501, LM_WRAPPER,FT_DEFAULT,0,LLL_VERBOSE);
	time.stop();
  

	os << "   fplll: " << time << endl << endl ;
   
	transpose(A,AT);
	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T12(A,NO_TRANSFORM,DEF_REDUCTION); //* name
	T12.isreduced(delta-0.1); //* name
      } 
   
    }// End on runs, k loop


    // END 
    fb.close();


 
  return 0;
}
