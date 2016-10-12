/* arithmetic plugins test file  

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

#include "hlll.h"

using namespace hplll; 

int main(int argc, char *argv[])  {
 
  int status=0;
  
  int d; 
  double delta; 

  int difference; 

  int succeed=0;
  int nbtest=0;

  filebuf fb;
  iostream os(&fb);

  //  *****************************************************  
  cout <<  "Testing arithmetic plug-ins" << endl; 
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
  cout << "     mpfr test " << endl; 

  d=4; 
  delta=0.99;

  A.resize(d+1,d); 
  AT.resize(d,d+1);  
  AT.gen_intrel(6);
  transpose(A,AT);

  
  Lattice<integer_t, mpfr_t, MatrixZT, matrix<FP_NR<mpfr_t> > > B1(A,NO_TRANSFORM,DEF_REDUCTION);
  B1.hlll(delta);

  transpose(AT,B1.getbase());

  TT.resize(d,d+1);
  fb.open ("416_in",ios::in);
  os >> TT ;
  fb.close();
  
  difference = !matcmp(AT, TT, d+1, d);

  status |= difference;
  
  if (difference) {
    cerr << "*** Invalid matrix comparison in arithmetic test" << endl;
  }
  else 
    succeed+=1;

  //  -------------------- TEST i --------------------------------

  nbtest+=1;
  cout << "     double test " << endl; 

  d=4; 
  delta=0.99;

  A.resize(d+1,d); 
  AT.resize(d,d+1);  
  AT.gen_intrel(6);
  transpose(A,AT);

  Lattice<integer_t, double, MatrixZT, matrix<FP_NR<double> > > B2(A,NO_TRANSFORM,DEF_REDUCTION);
  B2.hlll(delta);
  transpose(AT,B2.getbase());
  
  TT.resize(d,d+1);
  fb.open ("2_in",ios::in);
  os >> TT ;
  fb.close();

  difference = !matcmp(AT, TT, d+1, d);

  status |= difference;
  
  if (difference) {
    cerr << "*** Invalid matrix comparison in arithmetic test" << endl;
  }
  else 
    succeed+=1;

  //  -------------------- TEST i --------------------------------
#ifdef HPLLL_WITH_LONG_DOUBLE
  nbtest+=1;
  cout << "     long double test" << endl; 

  d=4; 
  delta=0.99;

  A.resize(d+1,d); 
  AT.resize(d,d+1);  
  AT.gen_intrel(6);
  transpose(A,AT);

  Lattice<integer_t, long double,  MatrixZT, matrix<FP_NR<long double> > > B3(A,NO_TRANSFORM,DEF_REDUCTION);
  B3.hlll(delta);
  transpose(AT,B3.getbase());
  TT.resize(d,d+1);
  fb.open ("3_in",ios::in);
  os >> TT ;
  fb.close();

  difference = !matcmp(AT, TT, d+1, d);

  status |= difference;
  
  if (difference) {
    cerr << "*** Invalid matrix comparison in arithmetic test" << endl;
  }
  else 
    succeed+=1;
#endif 


  //  *****************************************************  
  cout << endl << "     " << succeed << " arithmetic tests ok over " << nbtest << endl << endl; 
  //  *****************************************************  

  return status;
}
