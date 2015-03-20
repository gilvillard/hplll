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

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {

    // Reading of the zero 

  filebuf fb;
  iostream os(&fb);

  Z_NR<mpz_t> zz;
  int deg;

  fb.open ("zero120",ios::in);
  //  os >> bits ;
  //os >> n;
  os >> deg;
  os >> zz;
  fb.close();

  
  int setprec=34000;
  mpfr_set_default_prec(setprec);

  FP_NR<mpfr_t> zzf,tf;
  tf=10.0;
  set_z(zzf,zz);

  while (zzf.cmp(10.0) > 0) {
    zzf.div(zzf,tf);
  } 

  fb.open ("zero",ios::out);
  //  os >> bits ;
  //os >> n;
  os << zzf; 
  fb.close();
  
 
  // Real matrix 
  // -----------

 
  matrix<FP_NR<mpfr_t> > F;

  F.resize(1,deg+1);

  F.set(0,0,1.0);
  F.set(0,1,zzf);

  for(int k=2; k<deg+1; k++) {
    tf.mul(F.get(0,k-1),F.get(0,1));
    F.set(0,k,tf);
  }

 
  
  ZZ_mat<mpz_t> C;
 
  int start=utime();

  // Alpha must be less than prec by a factor of ||F|| for having alpha bits 
  relation_lift<long, double>(C, F, 30400, 800, FPLLL);
  
  //nullity=relations_lll<mpz_t, dpe_t, MatrixPE<double, dpe_t> > (C, F, setprec, 29200, 0);

  start = utime()-start;

  //print2maple(C,6,1);

  //cout << endl << "   Nullity: " << nullity << endl;
  cout << endl << "   Time: " << start/1000 << " ms" << endl;
 


  return 0;
}
