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



#ifdef _OPENMP
#include <omp.h>
#endif

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A0,A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------

  int n,d;
  double delta;
  
  command_line_basis(A0, n, d, delta, argc, argv); 
  
  A.resize(n,d);
  AT.resize(d,n);
  transpose(AT,A);

  
  OMPTimer time;

  {
    time.start();
  
#pragma omp parallel num_threads(2)
    {
      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A0,NO_TRANSFORM);
      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A0,NO_TRANSFORM);
      
      if (omp_get_thread_num() ==0)
	B.hlll(delta);
      
      if (omp_get_thread_num() ==1)
	C.hlll(delta);
      
    }
    
    time.stop();
    
    
    cout << "-- Parallel time:" << time  << endl; 
  }
    
     
   
   {
    time.start();

    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A0,NO_TRANSFORM);
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A0,NO_TRANSFORM);
      
    
    B.hlll(delta);
      
    
    C.hlll(delta);
      
    
    time.stop();
    
    
    cout << "-- Seq time:" << time  << endl; 
   }

   

  return 0;
}
