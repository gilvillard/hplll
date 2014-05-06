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
#include "tools.h" 

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------
  { 
  
    char ltype[1];
    char knapsack[]="r";
    char ajtai[]="a";
    char ntru[]="n";

    int d=8;
    int n;
    double delta = 0.75;
    int nbbits=10;
    double alpha=1.4;
    int q;

    PARSE_MAIN_ARGS {
      MATCH_MAIN_ARGID("-ltype",ltype);
      MATCH_MAIN_ARGID("-d",d);
      MATCH_MAIN_ARGID("-delta",delta);
      MATCH_MAIN_ARGID("-bits",nbbits);
      MATCH_MAIN_ARGID("-alpha",alpha);
      MATCH_MAIN_ARGID("-q",q);
      SYNTAX();
    }

    // Knapsack 
    // --------
    if (strcmp(ltype,knapsack) ==0) {
      n=d+1;
      A.resize(d+1,d); 
      AT.resize(d,d+1);  
      AT.gen_intrel(nbbits);
      transpose(A,AT);
    } 

    // Ajtai 
    // -----
    if (strcmp(ltype,ajtai) ==0) {
      n=d;
      A.resize(d,d); 
      AT.resize(d,d);  
      AT.gen_ajtai(alpha);
      transpose(A,AT);
    } 

    // NTRU like 
    // ---------
    if (strcmp(ltype,ntru) ==0) {
      d=2*d;
      n=d;
      A.resize(d,d); 
      AT.resize(d,d);  
      AT.gen_ntrulike(nbbits,q);
      transpose(A,AT);
    } 


    d=300;
    n=300;

    A.resize(d,d); 
    AT.resize(d,d); 

    filebuf fb;
    fb.open ("challenge-300",ios::in);
    iostream os(&fb);
    os >> AT ;
    fb.close();
    transpose(A,AT);

    

    print2maple(A,n,d);

    int start,startsec;

    cout << "--------------  HLLL" << endl << endl; 
    start=utime();
    startsec=utimesec();
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);
    B.hlll(delta);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << nbbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time A: " << start/1000 << " ms" << endl;
    cout << "   time A: " << startsec << " s" << endl;

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T1.isreduced(delta-0.1);

    cout << endl; 

    cout << "--------------  FPLLL" << endl << endl; 
    transpose(AT,A);

    start=utime();
    startsec=utimesec();
    lllReduction(AT, delta, 0.501, LM_FAST,FT_DEFAULT,0);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << nbbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time B: " << start/1000 << " ms" << endl;
    cout << "   time B: " << startsec << " s" << endl;

    transpose(A,AT);
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(A,NO_TRANSFORM,DEF_REDUCTION);
    T2.isreduced(delta-0.1);

    //print2maple(T1.getbase(),n,n);

    /*
    cout << endl; 
  
    transpose(AT,A);
 
    start=utime();
    startsec=utimesec();
    lllReduction(AT, delta, 0.5, LM_FAST,FT_DEFAULT,0);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << nbbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time C: " << start/1000 << " ms" << endl;
    cout << "   time C: " << startsec << " s" << endl;
  
    transpose(tmpmat,AT);
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T3(tmpmat,NO_TRANSFORM,DEF_REDUCTION);
    T3.isreduced(delta);
    */

  } 

 
  return 0;
}
