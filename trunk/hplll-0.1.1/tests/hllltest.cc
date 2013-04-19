/* hlll test file  

Created Dim 31 mar 2013 13:54:36 CEST 
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


#include "fplll.h"


using namespace std;
using namespace fplll;

#include "hlll.cc"  

int main(int argc, char *argv[])  {
 
  
  int d; 
  double delta; 

  int difference; 

  int succeed=0;
  int nbtest=0;

  int nbbits;


  //  *****************************************************  
  cout <<  "Testing hlll" << endl; 
  //  *****************************************************  

  // gmp integers 
  // ************

  typedef mpz_t integer_t;
  typedef matrix<Z_NR<integer_t> > MatrixZT;

  ZZ_mat<integer_t> A; // For hpLLL 
  ZZ_mat<integer_t> AT;  // fpLLL  
  ZZ_mat<integer_t> TT;  // test file  

  //  -------------------- TEST i --------------------------------
 
  nbtest+=1;

  d=100;
  nbbits=80;
  cout << "     hlll test, dim = " << d <<", " << nbbits << " bits, ldpe vs exp-dpe " << endl; 

  delta=0.99;

  A.resize(d+1,d); 
  AT.resize(d,d+1);  
  AT.gen_intrel(nbbits);
  transpose(A,AT);

  
  Lattice<integer_t, ldpe_t, MatrixZT, matrix<FP_NR<ldpe_t> > > B1(A,NO_TRANSFORM,DEF_REDUCTION);
  B1.hlll(delta);
  
  transpose(AT,B1.getbase());

  Lattice<integer_t, dpe_t, MatrixZT, matrixexp<double, dpe_t> >  B2(B1.getbase(),NO_TRANSFORM,DEF_REDUCTION);
  B2.hlll(delta);

  cout << "     " << B1.nblov << " vs " << B2.nblov << " Lovasz tests" << endl; 
  TT.resize(d,d+1);
  transpose(TT,B2.getbase());

  // Output in a file 
  /*
  transpose(AT,B1.getbase());  
  //fb.open ("old_out",ios::out); // Howto ? Pb with gmp ??? 
  //TT >> os;
  //fb.close();
  outfile.open ("foo",ios_base::trunc);
  outfile << AT  << endl;  
  outfile.close(); */
  
  difference = !matcmp(AT, TT, d+1, d);
  if (difference) {
    cerr << "*** Invalid matrix comparison in hlll test" << endl;
  }
  else 
    succeed+=1;

  //  -------------------- TEST i --------------------------------
  nbtest+=1;

  d=180;
  nbbits=240;
  cout << "     hlll test, dim = " << d <<", " << nbbits << " bits,  double vs exp-ldpe " << endl; 

  delta=0.99;

  A.resize(d+1,d); 
  AT.resize(d,d+1);  
  AT.gen_intrel(nbbits);
  transpose(A,AT);

  
  Lattice<integer_t, double, MatrixZT, matrix<FP_NR<double> > > B3(A,NO_TRANSFORM,DEF_REDUCTION);
  B3.hlll(delta);
  
  transpose(AT,B3.getbase());

  Lattice<integer_t, ldpe_t, MatrixZT, matrixexp<long double, ldpe_t> >  B4(B3.getbase(),NO_TRANSFORM,DEF_REDUCTION);
  B4.hlll(delta);

  cout << "     " << B3.nblov << " vs " << B4.nblov << " Lovasz tests" << endl; 
  TT.resize(d,d+1);
  transpose(TT,B4.getbase());
  
  difference = !matcmp(AT, TT, d+1, d);
  if (difference) {
    cerr << "*** Invalid matrix comparison in hlll test" << endl;
  }
  else 
    succeed+=1;

//  -------------------- TEST i --------------------------------
  nbtest+=1;

  d=12;
  nbbits=2000;
  cout << "     hlll test, dim = " << d <<", " << nbbits << " bits, mpfr vs exp-dpe " << endl; 

  delta=0.99;

  A.resize(d+1,d); 
  AT.resize(d,d+1);  
  AT.gen_intrel(nbbits);
  transpose(A,AT);

  
  Lattice<integer_t, mpfr_t, MatrixZT, matrix<FP_NR<mpfr_t> > > B5(A,NO_TRANSFORM,DEF_REDUCTION);
  B5.hlll(delta);
  
  transpose(AT,B5.getbase());

  Lattice<integer_t, dpe_t, MatrixZT, matrixexp<double, dpe_t> > B6(B5.getbase(),NO_TRANSFORM,DEF_REDUCTION);
 
  B6.hlll(delta);

  cout << "     " << B5.nblov << " vs " << B6.nblov << " Lovasz tests" << endl; 
  TT.resize(d,d+1);
  transpose(TT,B6.getbase());
  
  difference = !matcmp(AT, TT, d+1, d);
  if (difference) {
    cerr << "*** Invalid matrix comparison in hlll test" << endl;
  }
  else 
    succeed+=1;

  //  *****************************************************  
  cout << "     " <<  succeed << " hlll tests ok over " << nbtest << endl; 
  //  *****************************************************  


  return 0;
}
