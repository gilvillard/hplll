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


#include "matgen.h"
//#include "relations.cc" 

#include "nullspace.h" 

using namespace hplll;

/* ***********************************************

          MAIN   

   ********************************************** */


int main(int argc, char *argv[])  {

  
  int d,n,bitsize;
  int start,startsec,tps=0,tps_sec=0;

  d=1; 
  n=4;
  bitsize=4000;


  ZZ_mat<mpz_t> A;
  A.resize(d,n);
  A.gen_uniform(bitsize);

  ZZ_mat<mpz_t> C;

  start = utime(); 
  startsec = utimesec(); 
  
  nullspace_lll<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > (C, A, 0, 10000, HLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;
  
  cout << "--- cputime nullspace is " << tps/1000 << " ms" << endl;
  cout << "                         " << tps_sec << " s" << endl;

  cout <<  "Nullspace bit size = " << maxbitsize(C) << endl; 
  
  return 0;
}