/* 128 bit integers test file  

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

using namespace hplll; 

int main(int argc, char *argv[])  {
 
  int status=0;


  matrix<Z_NR<__int128_t> > A;
   
  A.resize(2,2);

  A(0,0)=4;
  A(0,1)=-4;
  A(1,0)=12321;
  A(1,1)=-8765432;

 

  //  *****************************************************  
  cout <<  "Testing 128 bits integers " << endl; 
  //  *****************************************************  

  matprod_in(A,A);
  matprod_in(A,A);
  
  

  Matrix<FP_NR<double> > B;
  B.resize(2,2);

  B(0,0).set_z(A(0,0));
  B(0,1).set_z(A(0,1));
  B(1,0).set_z(A(1,0));
  B(1,1).set_z(A(1,1));



  Matrix<Z_NR<mpz_t> > C;
  C.resize(2,2);

  C(0,0).set_f(B(0,0));
  C(0,1).set_f(B(0,1));
  C(1,0).set_f(B(1,0));
  C(1,1).set_f(B(1,1));

  filebuf fb;
  iostream os(&fb);

  Matrix<Z_NR<mpz_t> > CC;
  CC.resize(2,2);
  
  fb.open ("r128int",ios::in);
  os >> CC ;
  fb.close();

  
  int  difference = !matcmp(C, CC, 2, 2);

  status |= difference;

  if (difference) 
    cerr << "*** Invalid matrix comparison in 128 bit test" << endl;


  cout << "status " << status << endl;
  
  return status;
}
