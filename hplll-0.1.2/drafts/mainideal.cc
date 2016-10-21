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

#include "ideal.h"


/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 



int main(int argc, char *argv[])  {
  
  mat_ZZ B;
  int d,dbits; 

  PARSE_MAIN_ARGS {
    MATCH_MAIN_ARGID("--dim",d);
    MATCH_MAIN_ARGID("--dbits",dbits);
    SYNTAX();
  }

  d=generate_svp(B,d,dbits);

  cout << "Dimension " << d << endl; 

  typedef FP_NR<mpfr_t>   RT;
  typedef Z_NR<mpz_t>  ZT;
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  filebuf fb;
  fb.open ("ntl.txt",ios::out);
  iostream os(&fb);
  os << B ;
  fb.close();


  fb.open ("ntl.txt",ios::in);
  os >> AT ;
  fb.close();


  // ---------------------------------------------------------------------
 
  { 
  
    cout << "************************************************************************** " << endl; 
   
    int start,startsec;

    double delta=0.8;

    A.resize(d,d); 
    tmpmat.resize(d,d); 
    
    transpose(A,AT);

    cout << "--------------  HLLL" << endl << endl; 
    start=utime();
    startsec=utimesec();
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);
    B.hlll(delta);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << d*dbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time A: " << start/1000 << " ms" << endl;
    cout << "   time A: " << startsec << " s" << endl;

    //Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    //T1.isreduced(delta);

    cout << endl; 

    cout << "--------------  FPLLL" << endl << endl; 
    transpose(AT,A);

    start=utime();
    startsec=utimesec();
    lllReduction(AT, delta, 0.5, LM_WRAPPER,FT_DEFAULT,0);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << d*dbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time B: " << start/1000 << " ms" << endl;
    cout << "   time B: " << startsec << " s" << endl;

    //transpose(tmpmat,AT);
    //Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(tmpmat,NO_TRANSFORM,DEF_REDUCTION);
    //T2.isreduced(delta);

  } 

 
  return 0;
}