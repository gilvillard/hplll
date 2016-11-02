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
 
  int n;
  
  int found;

  int difference; 

  int succeed=0;
  int nbtest=0;

  filebuf fb;
  iostream os(&fb);

  int nbrel;

  long setprec;

  //  *****************************************************  
  cout <<  "Relation finder" << endl; 
  //  *****************************************************  

  matrix<FP_NR<mpfr_t> > A;   // Input matrix 

  typedef mpz_t integer_t;
  
  ZZ_mat<integer_t> C;  // Output relations 
  ZZ_mat<integer_t> Ccheck;  // Output relations 

 
  //  -------------------- TEST i --------------------------------
  nbtest+=1;


  // matrix<FP_NR<mpfr_t> > A;   // Input matrix
  // ZZ_mat<mpz_t> C;

  // int r=8;
  // int s=8;
  // int n=r*s+1;

  // int setprec=2800;
  // mpfr_set_default_prec(setprec);

  // gen3r2s(A,n,r,s);


  

  ZZ_mat<mpz_t> AZ;

  fb.open ("alpha.in",ios::in);
  os >> setprec ;
  os >> n;
  //AZ.resize(1,n);
  os >> AZ;
  fb.close();


 static string s;
  
  filebuf fb;
  iostream os(&fb);

 
  FP_NR<mpfr_t> f;   // Input matrix 
  
 
  fb.open ("t.in",ios::in);

  os >> s;

  fb.close();

  cout << s << endl;

  mpfr_set_str (f.get_data(), s.c_str(), 10, GMP_RNDN);

  cout << f << endl; :w

  
  mpfr_set_default_prec(setprec);
 

  FP_NR<mpfr_t> tmp;
  A.resize(1,n);
  for (int j=0; j<n; j++) {
    set_z(tmp,AZ(0,j));
    tmp.mul_2si(tmp,-setprec);
    A.set(0,j,tmp);
  }


  print2maple(A,1,n);
  
  nbrel=1;
  cout << endl;

  //print2maple(A,1,n);

  verboseDepth=0;

  Timer time;

  // ---------------------
  time.start();
  
  //found = relation_f<long, double>(C, A,240,60,800,20);

  found = relation_lll<dpe_t, MatrixPE<double, dpe_t> >(C, A,setprec,80,20,FPLLL);
  

  time.stop();
  // ---------------------
  
  print2maple(C,n,1);

  cout << "Time : " << time << endl; 
  
  

  return 0;
}
