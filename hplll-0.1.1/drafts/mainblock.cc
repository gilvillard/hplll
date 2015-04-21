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

  int n,d,k,i,j,l;
  double delta;

  command_line_basis(A0, n, d, delta, argc, argv); 

  int nblocks=2;
  int pcent = 5;   // pc * d > 100 
  
 
 
  vector< ZZ_mat<mpz_t> > TA(nblocks);

  vector<int> dimb(nblocks);

  // Slices of the orginal problem
  // -----------------------------
  
  if ((d%nblocks) ==0)
    for (k=0; k<nblocks; k++)
      dimb[k]=d/nblocks;
  else {
    for (k=0; k<d%nblocks; k++)
      dimb[k]=d/nblocks+1;
    for (k=d%nblocks; k<nblocks; k++)
      dimb[k]=d/nblocks;
  }

  // Plus overlap, but in the first block  
  // ------------------------------------
  
  int recouv=d*pcent/100;

  TA[0].resize(n,dimb[0]);
  
  for (k=1; k<nblocks; k++)
    TA[k].resize(n,dimb[k]+recouv);

  
  cout << endl << dimb << "    " << recouv << endl; 
    
  // Assign the non overlapping part
  // -------------------------------
  
  j=0;
  for (k=0; k<nblocks; k++) {
    for (l=0; l<dimb[k]; l++) {
      for (i=0; i<n;i++)
	(TA[k])(i,l)=A0(i,j);
      j++;
    }
  }

  // Assign the overlapping part, but for the 1rst block 
  // ---------------------------------------------------
 
  // k=0;
  // for (l=0; l<recouv; l++) {
  //   for (i=0; i<n;i++)
  //     (TA[k])(i,dimb[k]+l)=A0(i,dimb[k]+l);
  // }
  

  // The first vectors to the other blocks 
  for (k=1; k<nblocks; k++) 
    for (l=0; l<recouv; l++) 
      for (i=0; i<n;i++)
	(TA[k])(i,dimb[k]+l)=A0(i,l);
    

  // Dimension updates
  // -----------------

  int newd =0;
  newd+=dimb[0];
  
  // But the first block 
  for (k=1; k<nblocks; k++) {
      dimb[k]+=recouv;
      newd+=dimb[k];
  }


  // for (k=0; k<nblocks; k++)
  //   print2maple(TA[k],n, dimb[k]);
  
  // Reduction of the slices
  // -----------------------
 
  
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
     
  cout << endl << "Parallel time: " << time << endl;


  // The generating set 
  // ------------------
  
  ZZ_mat<mpz_t> L;
  L.resize(n,newd);

  
  // Reconstruction of a generating set from reduced slices  
  // ------------------------------------------------------

  j=0;
  for (k=0; k<nblocks; k++) {
    for (l=0; l<dimb[k]; l++) {
      for (i=0; i<n;i++)
	L(i,j)= (TA[k])(i,l); 
      j++;
    }
  }

 
  
  // int block_ind=0;
  // j=0;

  // for (l=0; l<dimb[0]; l++) {
    
  //   for (k=0; k<nblocks; k++) 
  //     if (block_ind < dimb[k]) {

  // 	for (i=0; i<n; i++)
  // 	  L(i,j)= (TA[k])(i,block_ind); 
  // 	j++;
  //     }

  //   block_ind+=1;
  // }
    

     
     // Reduction of the generating set 
    // -------------------------------


    time.start();
    GLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > W(L, d);
   
    W.hlll(delta);

    time.stop();
    cout << endl << "Remaining time: " << time << endl;
    
    // //print2maple(W.getbase(),n,newd);   
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(W.getbase(),NO_TRANSFORM,NO_LONG);
    T1.isreduced(delta-0.1);   

      // ***************

      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A0,NO_TRANSFORM);

      time.start();
      C.hlll(delta);
      time.stop();

      cout << endl << endl << "Conventional: " << time << endl;

      Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(C.getbase(),NO_TRANSFORM,NO_LONG);
      T2.isreduced(delta-0.1);   


    //   //print2maple(W.getbase(),n,d);
      
  return 0;
}
