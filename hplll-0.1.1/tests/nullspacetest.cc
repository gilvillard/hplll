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


#include "nullspace.h" 

using namespace hplll;

/* ***********************************************

          MAIN   

   ********************************************** */


int main(int argc, char *argv[])  {

 
  int d,n;


  int difference; 

  int succeed=0;
  int nbtest=0;

  filebuf fb;
  iostream os(&fb);

  int bitsize; 
  int prec;


  //  *****************************************************  
  cout <<  "Testing integer nullspace" << endl; 
  //  *****************************************************  

  typedef mpz_t integer_t;
 
  ZZ_mat<integer_t> A;
  ZZ_mat<integer_t> C;

  ZZ_mat<integer_t> T;

//  -------------------- TEST i --------------------------------
 
  nbtest+=1;

  d=20; 
  n=60;
  bitsize=20;

  //d=4;
  //n=8;
  //bitsize=8;

  A.resize(d,n);
  A.gen_uniform(bitsize);

  cout << "     direct integer decomp, " << d <<" x " << n <<", " << bitsize << " bits, mpz exp-dpe " << endl; 

  nullspace_direct_decomp <mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >(C,A);
 
  T.resize(n,n-d);

  
  fb.open ("C4_out",ios::in);
  os >> T ;
  fb.close();

  //print2maple(A,d,n);  
  //print2maple(C,n,n-d);

  difference = !matcmp(T, C, n-d, n);
  if (difference) {
    cerr << "*** Invalid matrix comparison in hlll test" << endl;
  }
  else 
    succeed+=1;
  
  //  -------------------- TEST i --------------------------------
 
  nbtest+=1;

  d=4; 
  n=20;
  bitsize=300;

  A.resize(d,n);
  A.gen_uniform(bitsize);

  cout << "     direct integer decomp, " << d <<" x " << n <<", " << bitsize << " bits, mpz double " << endl; 


  nullspace_direct_decomp <mpz_t, double , matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > >(C,A);
 
 
  T.resize(n,n-d);

  fb.open ("C5_out",ios::in);
  os >> T ;
  fb.close();
  //cout << C << endl; 
  difference = !matcmp(T, C, n-d, n);
  if (difference) {
    cerr << "*** Invalid matrix comparison in hlll test" << endl;
  }
  else 
    succeed+=1;

  //  -------------------- TEST i --------------------------------
  
  nbtest+=1;
  
  d=2; 
  n=10;
  bitsize=1000;

  prec=8000;

  A.resize(d,n);
  A.gen_uniform(bitsize);

  cout << "     indirect integer decomp, " << d <<" x " << n <<", " << bitsize << " bits, mpz double " << endl; 

  nullspace_indirect_decomp<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double > > > (C, A, prec); 

  T.resize(n,n-d);

  fb.open ("C6_out",ios::in);
  os >> T ;
  fb.close();
  //cout << C << endl; 
  difference = !matcmp(T, C, n-d, n);
  if (difference) {
    cerr << "*** Invalid matrix comparison in hlll test" << endl;
  }
  else 
    succeed+=1; 


 //  -------------------- TEST i --------------------------------
  
  nbtest+=1;
  
  d=5; 
  n=20;
  bitsize=40;
  

  prec=100;

  A.resize(d,n);
  A.gen_uniform(bitsize);

  cout << "     hjls integer decomp, " << d <<" x " << n <<", " << bitsize << " bits, mpz mpfr" << endl; 

  nullspace_hjls<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > (C, A, prec); 
  
  T.resize(n,n-d);

  fb.open ("C7_out",ios::in);
  os >> T ;
  fb.close();
 
  difference = !matcmp(T, C, n-d, n);
  if (difference) {
    cerr << "*** Invalid matrix comparison in hlll test" << endl;
  }
  else 
    succeed+=1;

 
 //  *****************************************************  
  cout << "     " <<  succeed << " nullspace tests ok over " << nbtest << endl; 
  //  *****************************************************  


  return 0;
}
