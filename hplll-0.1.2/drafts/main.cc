/* integer relations test file  

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

#include "matgen.h"
#include "relations.h" 

using namespace hplll;

/* ***********************************************

          MAIN   

   ********************************************** */


int main(int argc, char *argv[])  {
 
  int n;
  
  int found;

  int difference; 

  int succeed=0;
  int nbtest=0;

  int nbrel;

  long setprec;

  //  *****************************************************  
  cout <<  "Relation finder" << endl; 
  //  *****************************************************  

  matrix<FP_NR<mpfr_t> > A;   // Input matrix 

  typedef mpz_t integer_t;
  
  ZZ_mat<integer_t> C;  // Output relations 
  ZZ_mat<integer_t> Ccheck;  // Output relations 

 
  //  -------------------- TEST i --------------------------------
  nbtest+=1;


  // matrix<FP_NR<mpfr_t> > A;   // Input matrix
  // ZZ_mat<mpz_t> C;

  // int r=8;
  // int s=8;
  // int n=r*s+1;

  // int setprec=2800;
  // mpfr_set_default_prec(setprec);

  // gen3r2s(A,n,r,s);

  
  filebuf fb;
  iostream os(&fb);
  
  static string s;
  
  fb.open ("alpha.in",ios::in);

  os >> setprec ;
  os >> n;

  mpfr_set_default_prec(setprec);
  A.resize(1,n);
  for (int i=0; i<n; i++) {
    os >> s;
    mpfr_set_str (A(0,i).get_data(), s.c_str(), 10, GMP_RNDN);
  }

  fb.close();

 

  print2maple(A,1,n);
  
  nbrel=1;
  cout << endl;

  //print2maple(A,1,n);

  verboseDepth=0;

  Timer time;

  // ---------------------
  time.start();
  
  found = relation_f<long,  long double>(C, A,setprec,60,100,10,HLLL); // pas long double avec fplll 

  //found = relation<__int128_t, dpe_t, MatrixPE<double, dpe_t> >(C, A,setprec,60,10,80, HLLL); // pas 128_t fplll 
  

  time.stop();
  // ---------------------
  
  print2maple(C,n,1);

  cout << "Time : " << time << endl; 
  
  

  return 0;
}
