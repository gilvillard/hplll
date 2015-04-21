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
#include "matgen.h"


/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll;

int main() {

  int d=300;
  int nbbits=20;
 
  ZZ_mat<mpz_t> A0;
  
  matrix<Z_NR<mpz_t> > A;

  A0.resize(d,d); 
  A.resize(d,d);  
  A0.gen_uniform(nbbits);

  A.set(A0); 
 
  MatrixPE<double, dpe_t>  Af; 
  Af.resize(d,d);

  MatrixPE<double, dpe_t>  Bf; 
  Bf.resize(d,d);
 
  int j;

  for (j=0; j<d; j++) {
    Af.setcol(j,A.getcol(j),0,d);
    Bf.setcol(j,A.getcol(j),0,d);
  }

  OMPTimer time;

  FP_NR<dpe_t> x,y,one;
  x = 22.0;
  y = 22.0;
  one = 1.0;

  int NK=100;
  
  time.start();
  for (int k=0; k<NK; k++) {
    x.add(x,one);
    Af.submulcol(10, 20, x, d);
    Bf.submulcol(10, 20, x, d);
  } 
  time.stop();

  cout << endl << "**** seq: " << time << endl; 

  
  time.start();
  
#pragma omp parallel for num_threads(2)  shared(NK,d)
  for (int i =0; i<2; i++) {
    
    if (i ==0)
      for (int k=0; k<NK; k++) {
	 x.add(x,one);
	Af.submulcol(10, 20, x, d);
	
      }
    else
      for (int k=0; k<NK; k++) {
	 y.add(y,one);
	Bf.submulcol(10, 20, y, d);
      }
  }
  time.stop();

  cout << endl << "**** par: " << time << endl; 


  
  
}


