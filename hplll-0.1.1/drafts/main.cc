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
#include "hlllg.h"
#include "matgen.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A0,A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------

  int n,d,k,i,j;
  double delta;

  command_line_basis(A0, n, d, delta, argc, argv); 

  int nblocks=2;
  int recouv=2;
  
  vector< ZZ_mat<mpz_t> > TA(nblocks);

  //print2maple(A0,n,d);
    
  if ((d%nblocks) ==0) {
    for (k=0; k<nblocks-1; k++) {
      
      (TA[k]).resize(n,d/nblocks+recouv);

      for (i=0; i<n; i++)
	for (j=0; j<d/nblocks+recouv; j++)
	  (TA[k])(i,j)=A0(i,k*(d/nblocks)+j);
      
    
    }
    (TA[k]).resize(n,d/nblocks);

    for (i=0; i<n; i++)
      for (j=0; j<d/nblocks; j++)
	(TA[k])(i,j)=A0(i,k*(d/nblocks)+j);
      
    
    
  }
  else {
    for (k=0; k<nblocks-1; k++) {
      
      (TA[k]).resize(n,d/nblocks+1+recouv);

      for (i=0; i<n; i++)
	for (j=0; j<d/nblocks+1+recouv; j++)
	  (TA[k])(i,j)=A0(i,k*(d/nblocks+1)+j);
      
    
    }
    (TA[k]).resize(n,d/nblocks);

    for (i=0; i<n; i++)
      for (j=0; j<d/nblocks; j++)
	(TA[k])(i,j)=A0(i,k*(d/nblocks+1)+j);
    
    
  } 

  int newd=0;
  for (k=0; k<nblocks; k++)
    newd+=(TA[k]).getCols();

 
  ZZ_mat<mpz_t> L;
  L.resize(n,newd);

  int kacc=0;

 
    

  OMPTimer time;

  time.start();
#pragma omp parallel num_threads(nblocks)
    {
      int id=omp_get_thread_num();
      
      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(TA[id],NO_TRANSFORM);
      
      B.hlll(delta);

      TA[id]=B.getbase();
	
    }
  time.stop();

    cout << "********* " << time << endl;


    kacc=0;
    
    for (k=0; k<nblocks; k++) {
      for (j=0; j<(TA[k]).getCols(); j++) {
	for (i=0; i<n; i++) {
	  L(i,kacc)=(TA[k])(i,j);
	}
	kacc+=1;
      }
    } 

    //print2maple(L,n,newd);
    // Reduction of the generating set 
    // -------------------------------


    time.start();
    GLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > W(L, d);

    

    W.hlll(delta);

    //print2maple(W.getbase(),n,newd);   
   

  return 0;
}
