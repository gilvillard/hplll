/* Integer matrix nullspace test file  

Created Lun 14 oct 2013 16:28:07 CEST  
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


/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
 

  int n=0,d=0; 
  int decal=1;
  int setprec=53;

  /* ----------------------------------------- */

  typedef mpz_t integer_t;

  /* --------------------------------- */

  ZZ_mat<integer_t> A; 
  ZZ_mat<integer_t> AT;  

  int output=0;

  if (strcmp(argv[decal],"-u")==0) {

    n=atoi(argv[decal+1]);
    d=n;

    A.resize(n,d); 
    AT.resize(d,n);  

    AT.gen_uniform(atoi(argv[decal+2]));

    transpose(A,AT);

    output = atoi(argv[decal+3]);
    }

  if (strcmp(argv[decal],"-k")==0) {

    n=atoi(argv[decal+1])+1;
    d=n-1;

    A.resize(n,d); 
    AT.resize(d,n);  

    AT.gen_intrel(atoi(argv[decal+2]));
    transpose(A,AT);
    output = atoi(argv[decal+3]);

    cout << " " << endl;
  }

  if (strcmp(argv[decal],"-kk")==0) {

    d=atoi(argv[decal+1]);
    n=d+d/2;

    A.resize(n,d); 
    AT.resize(d,n);  

    for (int i=0; i<d/2; i++) 
      for (int j=0; j<d; j++) 
	A(i,j).randb(atoi(argv[decal+2]));

    for (int i=0; i<d; i++) 
      for (int j=0; j<d; j++) 
	A(i+d/2,j)=0; 

    for (int i=0; i<d; i++) 
      A(i+d/2,i)=1; 

    //print2maple(A,n,d);

    transpose(AT,A);
    output = atoi(argv[decal+3]);

    cout << " " << endl;
  }


  if (strcmp(argv[decal],"-a")==0) {

    n=atoi(argv[decal+1]);
    d=n;

    A.resize(n,d); 
    AT.resize(d,n);  

    AT.gen_ajtai(atof(argv[decal+2]));
    transpose(A,AT);
    output = atoi(argv[decal+3]);
    
    cout << " " << endl;
  }

  if (argc >= 6) { 
    if (strcmp(argv[decal+4],"-prec")==0) {
  
      setprec = atoi(argv[decal+5]);

    }
  }
  /* ------------------------------------------------ */
  int startsec,start;

  // DELTA ***********************

  long double delta=0.99;

  // *****************************

  /* FPLLL */ 
  
  int tpsfplll;
  
  start = utime(); //cputime();
  startsec = utimesec(); //cputime();

  lllReduction(AT, delta, 0.5, LM_FAST,FT_DOUBLE,0,LLL_VERBOSE); 
  lllReduction(AT, delta, 0.5, LM_FAST,FT_LONG_DOUBLE,0,LLL_VERBOSE); 
  lllReduction(AT, delta, 0.5, LM_PROVED,FT_DPE,0,LLL_VERBOSE); 

  tpsfplll=utimesec()-startsec;
  start=utime()-start;

  cout << endl; 
  cout << "**** Phase fplll  " << endl;
  cout << "      cputime fplll  is " << start/1000 << " ms";
  cout << "      cputime fplll  is " << tpsfplll << " s" << endl;


  /* HLLL */
  /* ---- */

  // Phase dpe double  
  //-----------------
  int tps1;
  
  transpose(A,AT);  // Selon avec ou pas fplll 

  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);

  startsec = utimesec(); //cputime();
  start = utime(); //cputime();

  B.hlll(delta);  

  tps1=utimesec()-startsec;
  start=utime()-start;

  cout << endl; 
  cout << "**** Phase dpe double  " << endl;
  cout << "      #Lovasz tests = " << B.nblov << endl;
  cout << "      cputime hlll  is " << start/1000 << " ms";
  cout << "      cputime hlll  is " << tps1 << " s" << endl;


  // Phase dpe long double  
  //-----------------
  int tps2;
  
  Lattice<mpz_t, ldpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<long double, ldpe_t> > C(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);

 start = utime(); //cputime();
  startsec = utimesec(); //cputime();

  C.hlll(delta);  

  tps2=utimesec()-startsec;
  start=utime()-start;

  cout << endl; 
  cout << "**** Phase dpe long double  " << endl;
  cout << "      #Lovasz tests = " << C.nblov << endl; 
  cout << "      cputime hlll  is " << start/1000 << " ms";
  cout << "      cputime hlll  is " << tps2 << " s" << endl;
 

  // Phase mpfr   
  //-----------
  int tps3;

  setprec=80;
  
  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > D(C.getbase(),NO_TRANSFORM,DEF_REDUCTION);

  D.setprec(setprec);

 start = utime(); //cputime();
  startsec = utimesec(); //cputime();

  D.hlll(delta);  

start=utime()-start;
  tps3=utimesec()-startsec;

  cout << endl; 
  cout << "**** Phase mpfr with prec " << setprec << endl;
  cout << "      #Lovasz tests = " << D.nblov << endl; 
  cout << "      cputime hlll  is " << start/1000 << " ms";
  cout << "      cputime hlll  is " << tps3 << " s" << endl;

  // Phase mpfr   
  //-----------

  setprec=90;
  
  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > E(D.getbase(),NO_TRANSFORM,DEF_REDUCTION);

  E.setprec(setprec);
  
  start = utime(); //cputime();
  startsec = utimesec(); //cputime();
  
  E.hlll(delta);  
  
  start=utime()-start;
  tps3=utimesec()-startsec;
  
  cout << endl; 
  cout << "**** Phase mpfr with prec " << setprec << endl;
  cout << "      #Lovasz tests = " << D.nblov << endl; 
  cout << "      cputime hlll  is " << start/1000 << " ms";
  cout << "      cputime hlll  is " << tps3 << " s" << endl;



  // Check 
  //------

  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > F(E.getbase(),NO_TRANSFORM,DEF_REDUCTION);

  unsigned int tps_check=0;

  cout << endl <<  "**** Phase check  "  << endl;

  start = utime(); //cputime();
  startsec = utimesec();
  
  F.isreduced(0.8);
  
  start=utime()-start;
  tps_check=utimesec()-startsec;
  
  cout << endl; 
  
  cout << "      cputime hlll  is " << start/1000 << " ms" << endl;
  cout << "      cputime check  is " << tps_check << " s" << endl;

  if (output) print2maple(D.getbase(),n,d);

  return 0;
}
