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
#include "matgen.h"
#include "relations.h"

#include "tools.h"

using namespace hplll;


/* ***********************************************

          MAIN   

   ********************************************** */



int main(int argc, char *argv[])  {
  
  // TEST BRUNO
  // Types
  // *****

  typedef mpfr_t RT;
  typedef double FT;


  filebuf fb;
  iostream os(&fb);

  ZZ_mat<mpz_t> A;
  int bits, n;
  fb.open ("tbs.in",ios::in);
  os >> bits ;
  os >> n;
  A.resize(1,n);
  os >> A;
  fb.close();


  int setprec=bits;
  mpfr_set_default_prec(setprec);
  cout << "Bits: " << bits << ", " <<  n << " input real numbers" << endl << endl;


  FP_NR<RT> quodigits;
  quodigits=2;
  quodigits.pow_si(quodigits,-bits);

  quodigits.print();

  matrix<FP_NR<RT> > F;   // Input matrix
  FP_NR<RT> tmp;
  F.resize(1,n);
  for (int j=0; j<n; j++) {
    set_z(tmp,A(0,j));
    tmp.mul(tmp,quodigits);
    F.set(0,j,tmp);
  }

  print2maple(F,1,n);

  ZZ_mat<mpz_t> C;


  int start = utime();

  relation_lift<long, double>(C, F, setprec, 800, FPLLL);

  
  start = utime()-start;
  cout << endl << "   Time: " << start/1000 << " ms" << endl;
  cout << endl << "   Time: " << start << " us" << endl;


  return 0;
}
