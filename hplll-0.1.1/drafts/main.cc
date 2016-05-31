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


#include "hplll.h"

#include "slll.h"

//#include "slll-wrap.h"
#include "wrappers.h"

using namespace hplll; 


    
/* ***********************************************

          MAIN   

   ********************************************** */



int main(int argc, char *argv[])  {

  typedef mpz_t  ZT;
  //typedef long ZT;

  ZZ_mat<mpz_t> A0; // For hpLLL
  
  ZZ_mat<ZT> A; // For hpLLL 

  ZZ_mat<ZT> AT;  // fpLLL
 
  // ---------------------------------------------------------------------
   
  int d=8;
  int n;
  
  //int nbbits=8;
  //double alpha;
  //int output;
  
  double delta = 0.99;

  //int m=1;
  
  //int lovmax=1000000;

  int S=4;  // Rajouter à commandline 

  int nbthreads=4;

  command_line_basis(A0, n, d, delta, argc, argv);

  // Attention en 128 bits, mpfr get_si pas autrement 
  matrix_cast(A,A0);  // temporaire avec ci-dessous

  // Random unimodular preconditioning for the 512 example
  // -----------------------------------------------------

  ZZ_mat<ZT> Q;
  Q.resize(d,d);

  for (int i=0; i<d; i++)
    Q(i,i)=1;
  
  for (int i=0; i<d; i++)
    for (int j=i+1; j<d ; j++) 
      Q(i,j).randb(4);

  for (int i=0; i<d; i++)
    for (int j=0; j<i ; j++) 
      Q(i,j)=0;
 
  matprod_in(A,Q);
  
  for (int i=0; i<d; i++)
    Q(i,i)=1;
  
  for (int i=0; i<d; i++)
    for (int j=i+1; j<d ; j++) 
      Q(i,j)=0;
  
   for (int i=0; i<d; i++)
    for (int j=0; j<i ; j++) 
      Q(i,j).randb(4);

  matprod_in(A,Q);
  
  
  
  // With the wrapper
  // ----------------
  
  ZZ_mat<ZT> C; 

  Timer tw;
 
  verboseDepth=3;

  cout  << "Dimension " << n << "     " << d << endl; 

  tw=slll<mpz_t>(C,A,255,4,delta, true, true); 
    
  //tw.start();
  // slll_wrap<mpz_t, ldpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<long double, long dpe_t> > (C,A,255,4,delta,SEYSEN_REDUCTION);
  //tw.stop();
  
  //cout << endl << "lllw: " << tw << endl; // cf bout de hlll aussi pour l'instant 

  //cout << "Transposed result" << endl;

  //cout << transpose(C) << endl; 
  
   // With hlll
   // ---------

  verboseDepth=0;
   
  Lattice<ZT, ldpe_t, matrix<Z_NR<ZT> >, MatrixPE<long double, ldpe_t> > L(A,NO_TRANSFORM,DEF_REDUCTION);
   
   Timer tl;
   tl.start();
   //L.hlll(delta);
   tl.stop();

  
   cout << endl << "hlll: " << tl << endl;

   // Verification
   // -------------

   cout << endl << endl;
  
   Lattice<ZT, mpfr_t,  matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > Btest(L.getbase(),NO_TRANSFORM,DEF_REDUCTION);
   Btest.isreduced(delta-0.1);

   //Lattice<ZT, mpfr_t,  matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > Ltest(L.getbase(),NO_TRANSFORM,DEF_REDUCTION);
   //Ltest.isreduced(delta-0.1);


   //DBG ratio 
   double t,u,v,w;

   ratio<ZT>(C,t,u,v,w);
   
   cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
   cout << ".. Average diagonal ratio: " << u << endl;
   cout << ".. Max diagonal ratio: " << v << endl;
   cout << ".. First vector quality: " << w << endl;

   cout << "-----------------------" << endl;

   cout << "SLLL: " << tw << endl;
   //tw.print(cout);
   cout << "HLLL :" << tl << endl;
   //tl.print(cout);  

    
  return 0;
}
