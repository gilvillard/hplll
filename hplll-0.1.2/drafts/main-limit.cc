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
  int method=0;

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
    if (output==1)  print2maple(A,n,d);
  }

  if (strcmp(argv[decal],"-k")==0) {
    
    n=atoi(argv[decal+1])+1;
    d=n-1;

    A.resize(n,d); 
    AT.resize(d,n);  

    AT.gen_intrel(atoi(argv[decal+2]));
    transpose(A,AT);
    output = atoi(argv[decal+3]);

    if (output==1)  print2maple(A,n,d);
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

    transpose(AT,A);
    output = atoi(argv[decal+3]);

    if (output==1)  print2maple(A,n,d);

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

    if (output == 1) print2maple(A,n,n);
  }

  if (argc >= 5) {   
    method = atoi(argv[decal+4]);
    
  }
 
  /* ------------------------------------------------ */
  int startsec,start;

  // DELTA ***********************

  long double delta=0.99;

  // *****************************
  /* FPLLL WRAPPER */ 
  
  if (method == -2) {
    int tpsfplll;
  
    start = utime(); //cputime();
    startsec = utimesec(); //cputime();
    
    lllReduction(AT, delta, 0.5, LM_WRAPPER,FT_DEFAULT,0,LLL_VERBOSE); 

    tpsfplll=utimesec()-startsec;
    start=utime()-start;
    
    cout << endl; 
    cout << "**** Phase fplll  " << endl;
    cout << "        cputime fplll  is " << start/1000 << " ms";
    cout << "        cputime fplll  is " << tpsfplll << " s" << endl;

    transpose(A,AT);

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > E(A,NO_TRANSFORM,DEF_REDUCTION);

    unsigned int tps_check=0;
    
    cout << endl <<  "**** Phase check  "  << endl;
    
    start = utime(); //cputime();
    startsec = utimesec();
    
    E.isreduced(0.8);
    
    start=utime()-start;
    tps_check=utimesec()-startsec;
    
    cout << endl; 
    
    cout << "      cputime check  is " << start/1000 << " ms" << endl;
    cout << "      cputime check  is " << tps_check << " s" << endl;
  }

  /* FPLLL HEURISTIC */ 
  
  if (method == -1) {
    int tpsfplll;
  
    start = utime(); //cputime();
    startsec = utimesec(); //cputime();
    
    lllReduction(AT, delta, 0.5, LM_HEURISTIC,FT_DPE,0,LLL_VERBOSE); 
   

    tpsfplll=utimesec()-startsec;
    start=utime()-start;
    
    cout << endl; 
    cout << "**** Phase fplll  " << endl;
    cout << "        cputime fplll  is " << start/1000 << " ms";
    cout << "        cputime fplll  is " << tpsfplll << " s" << endl;

    transpose(A,AT);

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > E(A,NO_TRANSFORM,DEF_REDUCTION);

    unsigned int tps_check=0;
    
    cout << endl <<  "**** Phase check  "  << endl;
    
    start = utime(); //cputime();
    startsec = utimesec();
    
    E.isreduced(0.8);
    
    start=utime()-start;
    tps_check=utimesec()-startsec;
    
    cout << endl; 
    
    cout << "      cputime check  is " << start/1000 << " ms" << endl;
    cout << "      cputime check  is " << tps_check << " s" << endl;
  }

  /* FPLLL PROVED */ 
  
  if (method == 0) {
    int tpsfplll;
  
    start = utime(); //cputime();
    startsec = utimesec(); //cputime();
    
    lllReduction(AT, delta, 0.5, LM_PROVED,FT_DPE,0,LLL_VERBOSE); 
    //lllReduction(AT, delta, 0.5, LM_WRAPPER,FT_DPE,0,LLL_VERBOSE); 

    tpsfplll=utimesec()-startsec;
    start=utime()-start;
    
    cout << endl; 
    cout << "**** Phase fplll  " << endl;
    cout << "        cputime fplll  is " << start/1000 << " ms";
    cout << "        cputime fplll  is " << tpsfplll << " s" << endl;

    transpose(A,AT);

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > E(A,NO_TRANSFORM,DEF_REDUCTION);

    unsigned int tps_check=0;
    
    cout << endl <<  "**** Phase check  "  << endl;
    
    start = utime(); //cputime();
    startsec = utimesec();
    
    E.isreduced(0.8);
    
    start=utime()-start;
    tps_check=utimesec()-startsec;
    
    cout << endl; 
    
    cout << "      cputime check  is " << start/1000 << " ms" << endl;
    cout << "      cputime check  is " << tps_check << " s" << endl;
  }

  /* HLLL */
  /* ---- */
  if (method==1) {

    int tps1;
    
    transpose(A,AT);  // Selon avec ou pas fplll 
    
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);
    
    startsec = utimesec(); //cputime();
    start = utime(); //cputime();
    
    B.hlll(delta,LLL_VERBOSE);  
    
    tps1=utimesec()-startsec;
    start=utime()-start;
    
    cout << endl; 
    
    cout << "      #Lovasz tests = " << B.nblov << endl;
    cout << "      cputime hlll  is " << start/1000 << " ms";
    cout << "      cputime hlll  is " << tps1 << " s" << endl;


    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > E(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);

    unsigned int tps_check=0;

    cout << endl <<  "**** Phase check  "  << endl;
    
    start = utime(); //cputime();
    startsec = utimesec();
    
    E.isreduced(0.8);
  
    start=utime()-start;
    tps_check=utimesec()-startsec;
  
    cout << endl; 
  
    cout << "      cputime check  is " << start/1000 << " ms" << endl;
    cout << "      cputime check  is " << tps_check << " s" << endl;
    
   
  }
  /* SEYSEN */
  /* ------ */
  if (method==2) {

    int tps1;
    
    transpose(A,AT);  // Selon avec ou pas fplll 
    
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,SEYSEN_REDUCTION);
    
    startsec = utimesec(); //cputime();
    start = utime(); //cputime();
    
    B.hlll(delta);  
    
    tps1=utimesec()-startsec;
    start=utime()-start;
    
    cout << endl; 
    
    cout << "      #Lovasz tests = " << B.nblov << endl;
    cout << "      cputime seysen  is " << start/1000 << " ms";
    cout << "      cputime seysen  is " << tps1 << " s" << endl;


    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > E(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);

    unsigned int tps_check=0;

    cout << endl <<  "**** Phase check  "  << endl;
    
    start = utime(); //cputime();
    startsec = utimesec();
    
    E.isreduced(0.8);
  
    start=utime()-start;
    tps_check=utimesec()-startsec;
  
    cout << endl; 
  
    cout << "      cputime check  is " << start/1000 << " ms" << endl;
    cout << "      cputime check  is " << tps_check << " s" << endl;
    
   
  }

  return 0;
}
