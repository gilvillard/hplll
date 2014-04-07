/* integer relations test file  

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

#include "matgen.h"
#include "relations.h" 

using namespace hplll;

/* ***********************************************

          MAIN   

   ********************************************** */


int main(int argc, char *argv[])  {
 
  int n,d;
  int r,s,t,u;

  int difference; 

  int succeed=0;
  int nbtest=0;

  filebuf fb;
  iostream os(&fb);

  int nbrel;

  long setprec;

  //  *****************************************************  
  cout <<  "Testing relation finder" << endl; 
  //  *****************************************************  

  matrix<FP_NR<mpfr_t> > A;   // Input matrix 

  typedef mpz_t integer_t;
  typedef matrix<Z_NR<integer_t> > MatrixZT;

  ZZ_mat<integer_t> C;  // Output relations 
  ZZ_mat<integer_t> Ccheck;  // Output relations 

  //  -------------------- TEST i --------------------------------
  nbtest+=1;

  r=4; 
  s=4; 
  n=r*s+1;
  d=1;
  
  
  setprec=320;
  mpfr_set_default_prec(setprec);

  gen3r2s(A,n,r,s);

  nbrel=1;

  cout << "     hjls test matrix double, dim = " << n <<", " << setprec << " bits " << endl; 

  nbrel=relations_hjls<mpfr_t,double,matrix<FP_NR<double> > >(C,A,nbrel,setprec);

  Ccheck.resize(n,1);
  fb.open ("C1_out",ios::in);
  os >> Ccheck ;
  fb.close();

  if (nbrel==1) {
    difference = !matcmp(C, Ccheck, 1, n);
    if (difference) {
      cerr << "*** Invalid matrix comparison in hjls test" << endl;
    }
    else 
      succeed+=1;
  }

#ifdef HPLLL_WITH_LONG_DOUBLE
   //  -------------------- TEST i --------------------------------
  nbtest+=1;

  r=5; 
  s=5; 
  n=r*s+1;
  d=1;
  
  
  setprec=600;
  mpfr_set_default_prec(setprec);

  gen3r2s(A,n,r,s);

  nbrel=1;

  cout << "     hjls test exp long dpe, dim = " << n <<", " << setprec << " bits " << endl; 

  nbrel=relations_hjls<mpfr_t,ldpe_t, MatrixPE<long double, ldpe_t> >(C,A,nbrel,setprec);
 

  Ccheck.resize(n,1);
  fb.open ("C2_out",ios::in);
  os >> Ccheck ;
  fb.close();

  if (nbrel==1) {
    difference = !matcmp(C, Ccheck, 1, n);
    if (difference) {
      cerr << "*** Invalid matrix comparison in hjls test" << endl;
    }
    else 
      succeed+=1;
  }
#endif 
  
  //  -------------------- TEST i --------------------------------
  nbtest+=1;

  r=2;
  s=3;
  t=3;
  u=2;
  n=r*s+t*u+1;

  setprec=660;
  mpfr_set_default_prec(setprec);

  gen3r2s7t5u(A,n,r,s,t,u); 

  nbrel=1;

  cout << "     hjls 2 vector test exp dpe, dim = " << n <<", " << setprec << " bits " << endl; 

  nbrel=relations_hjls<mpfr_t, dpe_t, MatrixPE<double, dpe_t> >(C,A,nbrel,setprec);
 
  Ccheck.resize(n,1);
  fb.open ("C3_out",ios::in);
  os >> Ccheck ;
  fb.close();

  if (nbrel==1) {
    difference = !matcmp(C, Ccheck, 1, n);
    if (difference) {
      cerr << "*** Invalid matrix comparison in hjls test" << endl;
    }
    else 
      succeed+=1;
  }
  

  //  *****************************************************  
  cout << "     " << succeed << " relations tests ok over " << nbtest << endl; 
  //  *****************************************************  

  return 0;
}
