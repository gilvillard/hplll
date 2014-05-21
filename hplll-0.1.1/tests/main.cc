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
#include "plll.h"



/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT;  // fpLLL  

  // ---------------------------------------------------------------------
  
    cout << "************************************************************************** " << endl; 

    int n,d;
    double delta;
    
    
    command_line_basis(A, n, d, delta, argc, argv); 

      
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,TRANSFORM,DEF_REDUCTION);

    int i;

    for (i=0; i<n; i++) {
      B.householder_r(i);
      B.householder_v(i);
    }

    for (i=n; i<d; i++) {
      B.householder_r(i);
      B.hsizereduce(i,n-1);
     
      }
   
 
    print2maple(B.getR(),n,d);
    print2maple(B.getU(),d,d);

  return 0;
}
