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


#include "matgen.cc"
#include "relations.cc" 

#include "nullspace.cc" 


#include "matmixed.h" 

#include "matpeomp.h" 




/* ***********************************************

          MAIN   

   ********************************************** */


int main() {

  int d=800;
  int nbbits=200;
 
  ZZ_mat<mpz_t> A0; 
  matrix<Z_NR<mpz_t> > A;

  A0.resize(d,d); 
  A.resize(d,d);  
  A0.gen_uniform(nbbits);

  A.set(A0); 
 
  MatrixPE<double, dpe_t>  Af; 
  Af.resize(d,d);

  MatrixPE_omp<double, dpe_t>  Afomp; 
  Afomp.resize(d,d);

  int j;

  for (j=0; j<d; j++) {
    Af.setcol(j,A.getcol(j),0,d);
    Afomp.setcol(j,A.getcol(j),0,d);
  }

  //print2maple(A,d,d);


  int k;   

  FP_NR<dpe_t> a;
  a=10.0;
  
  int K=1;

  // Single 
  int start=utime();
  for (k=0; k<K; k++) 

    Af.submulcol(2,d-1,a,d);

  start = utime()-start;
  cout << endl << "   Time: " << start << " us" << endl;

  // OMP 
  start=utime();
  for (k=0; k<K; k++) 

    Afomp.submulcol(2,d-1,a,d);

  start = utime()-start;
  cout << endl << "  OMP Time: " << start << " us" << endl;

  //cout << endl << "   OMP Time: " << start/1000 << " ms" << endl;
//printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
 
}


