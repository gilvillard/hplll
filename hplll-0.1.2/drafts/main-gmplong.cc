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
#include "lehmer.cc"

#include "tools.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {


  int d=8;
  int nbbits = 10;
  
  PARSE_MAIN_ARGS {
    MATCH_MAIN_ARGID("-d",d);
    MATCH_MAIN_ARGID("-bits",nbbits);
    SYNTAX();
  }

  ZZ_mat<mpz_t> A,B,C;

  int i,j;
  
  A.resize(d,d);
  B.resize(d,d);
  C.resize(d,d);

  for (i=0; i<d; i++)
    for (j=0; j<d ; j++) {
      A(i,j).randb(nbbits);
      B(i,j).randb(nbbits);
    }


  matrix<Z_NR<long> > Ap,Bp,Cp;
  Ap.resize(d,d);
  Bp.resize(d,d);
  Cp.resize(d,d);

  for (i=0; i<d; i++)
  for (j=0; j<d ; j++) {
    Ap(i,j)=mpz_get_si(A(i,j).get_data());
    Bp(i,j)=mpz_get_si(B(i,j).get_data());
  }


  ZZ_mat<double> Af,Bf,Cf;
  Af.resize(d,d);
  Bf.resize(d,d);
  Cf.resize(d,d);

  for (i=0; i<d; i++)
  for (j=0; j<d ; j++) {
    Af(i,j)=((double) Ap(i,j).get_data());
    Bf(i,j)=((double) Bp(i,j).get_data());
  }

    
  int start;
  
  start=utime();
  
  matprod(C,A,B);

  start=utime()-start;

  cout << "    time mpz : " << start/1000 << " ms" << endl;

  start=utime();
  
  matprod(Cp,Ap,Bp);

  start=utime()-start;

  cout << "    time long : " << start/1000 << " ms" << endl;

  start=utime();
  
  matprod(Cf,Af,Bf);

  start=utime()-start;

  cout << "    time double : " << start/1000 << " ms" << endl;

  
  
  return 0;
}
