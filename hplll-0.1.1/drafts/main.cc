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

  ZZ_mat<mpz_t> AT;  // fpLLL
 
  // ---------------------------------------------------------------------
  { 
  
  int d=8;
  int n;
  
  int nbbits=8;
  double alpha;
  int output;
  
  double delta = 0.99;

  int K=8;

  int m=1;
  
  int lovmax=1000000;

  command_line_basis(A, n, d, delta, argc, argv);

  char type[]="";
    
  PARSE_MAIN_ARGS {
    MATCH_MAIN_ARGID("-type",type);
      MATCH_MAIN_ARGID("-d",d);
      MATCH_MAIN_ARGID("-m",m);
      MATCH_MAIN_ARGID("-delta",delta);
      MATCH_MAIN_ARGID("-bits",nbbits);
      MATCH_MAIN_ARGID("-alpha",alpha);
      MATCH_MAIN_ARGID("-lovmax",lovmax);
      MATCH_MAIN_ARGID("-output",output);
      MATCH_MAIN_ARGID("-K",K);
      SYNTAX();
    }

  
  int i,j;


  // Lecture de A partiellement réduite

  // filebuf fb;
  // iostream os(&fb);

  // //fb.open ("in78",ios::in);
  // //n=100;  K=78; d=79;
  
  // //fb.open ("in98",ios::in);
  // //n=100;  K=98; d=99;

  // fb.open ("in118",ios::in);
  // n=120;  K=118; d=119;

  
  // A.resize(n,d);
  // os >> A ;
  // fb.close();

  // Découpage rectangle de A
  
 
  ZZ_mat<mpz_t>  T;
  T.resize(n,K);
  
  for (i=0; i<n; i++)
    for (j=0; j<K; j++)
      T(i,j)=A(i,j);

  
  // Réduction du premier rectangle
  
  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > LT(T,NO_TRANSFORM,DEF_REDUCTION);
 
  Timer lt;
  lt.start();
  LT.hlll(delta);
  lt.stop();
  cout << endl << "d: " << d << "    " <<  "K: " << K << endl << endl; 
  cout << "1er bout: " << lt << "   " <<  LT.nbswaps << endl;
  
  // On complète de 1

  T=LT.getbase();
  
  ZZ_mat<mpz_t>  B;
  B.resize(n,K+1);
  
  for (i=0; i<n; i++)
    for (j=0; j<K; j++)
      B(i,j)=T(i,j);

   for (i=0; i<n; i++)
     B(i,j)=A(i,K);

 
   
   // Avec plll
   // ----------
   
   // PLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > LP(B,NO_TRANSFORM,DEF_REDUCTION);
   
   // Timer lp;
   // lp.start();
   // LP.hlll(delta);
   // lp.stop();
   // cout << endl << "plll: " << lp << "   " <<  LP.nbswaps << endl;
   // cout << endl << "compteur: " << LP.compteur <<  endl;

   Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > TT(B,NO_TRANSFORM,DEF_REDUCTION);

   Timer lp0,lp;

   

   lp0.start();
   TT.householder();
   TT.hsizereduce(K);
   lp0.stop();
   
  
   SLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t>  > LP(TT.getbase(),8,TRANSFORM,DEF_REDUCTION);
 
   lp.start();
   LP.hlll(delta,53,8,lovmax);
   lp.stop();
   cout << endl << "hsize: " << lp0  << endl;
   cout << endl << "slll: " << lp  << endl;

   // Avec hlll
   // ----------
   
   Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > LB(B,NO_TRANSFORM,DEF_REDUCTION);
   
   Timer lb;
   lb.start();
   LB.hlll(delta);
   lb.stop();

  
   cout << endl << "hplll: " << lb << "   " <<  LB.nbswaps << endl;
   cout << endl << "compteur: " << LB.compteur <<  endl;

   // Avec fplll
   // ----------

   ZZ_mat<mpz_t> BT;  
   BT.resize(K+1,n);
   transpose(BT,B);
   
   Timer fp;
   fp.start();
   lllReduction(BT, delta, 0.501, LM_WRAPPER);
   fp.stop();
   cout << endl << "fplll: " << fp << endl;

   
   cout << endl;
   
   Lattice<mpz_t, mpfr_t,  matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > Btest(LP.getbase(),NO_TRANSFORM,DEF_REDUCTION);
   Btest.isreduced(delta-0.1);

  } 

 
  return 0;
}
