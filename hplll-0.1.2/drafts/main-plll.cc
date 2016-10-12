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
#include "slll.h"

#include "tools.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  
  int n;
  int d=8;
  
  double delta = 0.99;

  command_line_basis(A, n, d, delta, argc, argv);

  //Timer lt;
  //lt.start();

  // Seysen usual
  // ------------
  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,SEYSEN_REDUCTION);
   
  Timer tusual;

  tusual.start();
  
  B.hlll(delta);

  tusual.stop();
  
  double t,u,v,w;
  ratio<mpz_t>(B.getbase(),t,u,v,w);

  cout << endl << ".. LLL Seysen, log 2 Frobenius norm cond: " << t << endl;
  cout << endl << "   Time usual: " << tusual << endl;

  // New Seysen
  // ----------
  
  cout << endl; 
  PLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > LP(A,NO_TRANSFORM,DEF_REDUCTION);

  Timer tnew;
 
  
  tnew.start();

  LP.hlll(0.99); 
    
  tnew.stop();
  
  ratio<mpz_t>(LP.getbase(),t,u,v,w);

  cout << endl << ".. New Seysen, log 2 Frobenius norm cond: " << t << endl;
  cout << endl << "   Time new: " << tnew << endl;

  // LP.householder(); 
  // LP.seysenreduce(0,d-1);
 
  // ratio<mpz_t>(LP.getbase(),t,u,v,w);

  // cout << endl << ".. New Seysen 2, log 2 Frobenius norm cond: " << t << endl << endl;

  //  LP.householder(); 
  // LP.seysenreduce(0,d-1);
 
  // ratio<mpz_t>(LP.getbase(),t,u,v,w);

  // cout << endl << ".. New Seysen 2, log 2 Frobenius norm cond: " << t << endl << endl;
   
  return 0;
}
