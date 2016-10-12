/* NTL tests   

Created Jeu  2 jui 2016 16:28:42 CEST  
Copyright (C) 2016      Gilles Villard 


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

#include "wrappers.h"


#include <NTL/LLL.h>


/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

using namespace NTL;

int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A;

  ZZ_mat<mpz_t> C;
 
 
  // ---------------------------------------------------------------------

  int n,d;
  double delta;

  command_line_basis(A, n, d, delta, argc, argv); 

 
  Timer th,tf,tn;
  
  // HLLL ------------------------------------------
   
  
  //Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);

  // verboseDepth = 1;
  // th.start();
  // status=B.hlll(delta);
  // th.stop();
  
  verboseDepth = 1;
  
  th=hlll<mpz_t>(C, A, 0.99, true, true);
    
  //th=hlll<__int128_t>(C, A, 0.99, true,true); 
  //hlll<long>(C, A, 0.99, false, true); 
 
  cout << endl; 

  filebuf fb;
  fb.open ("ntl.txt",ios::out);
  iostream os(&fb);
  os <<  transpose(A) ;
  fb.close();

  cout << "--------------  FPLLL " << endl << endl;


  ZZ_mat<mpz_t> AT;

  AT.resize(d,n);
  
  transpose(AT,A);

  tf.start();

  lllReduction(AT, delta, 0.501, LM_WRAPPER, FT_DEFAULT,0,LLL_VERBOSE);
 

  tf.stop();
  
  
  cout << "--------------  NTL " << endl << endl;


  Mat<ZZ> B; // (INIT_SIZE, 2,2) ;
 
  fb.open ("ntl.txt",ios::in);
  os >>  B ;
  fb.close();  
  system("rm ntl.txt");

  tn.start();

  LLL_XD(B,0.99); 
  
  tn.stop();

  cout << "----------------------- CHECK " << endl;
  
  fb.open ("ntlout.txt",ios::out);
  os <<  B ;
  fb.close();  

  ZZ_mat<mpz_t> TNC,NC;

  NC.resize(d,n);
  TNC.resize(n,d);

 
  fb.open ("ntlout.txt",ios::in);
  os >> NC ;
  fb.close();  
  system("rm ntlout.txt");

  transpose(TNC,NC);
  
  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(TNC,NO_TRANSFORM,DEF_REDUCTION);
  verboseDepth=0;
  T.isreduced(delta-0.1);

  //  cout << "-----------------------" << endl;

   cout << "HLLL: " << th << endl;

   cout << "FPLLL :" << tf << endl;
   
   cout << "NTL :" << tn << endl;
  

  return 0;
}
