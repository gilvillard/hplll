/* Some matrix constructors  

Created Sam  6 avr 2013 17:42:48 CEST 
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

#include  "matgen.h"

#ifndef MATGEN_CC
#define MATGEN_CC

namespace hplll {

/* ***********************************************

          BASIS GENERATION FROM THE COMMAND LINE    

   ********************************************** */

  void command_line_basis(ZZ_mat<mpz_t>& A, int& n, int& d, double &delta, int argc, char *argv[]) {

    ZZ_mat<mpz_t> AT;

    char type[]="r";
    char knapsack[]="r";
    char block[]="rb";
    char ajtai[]="a";
    char ntru[]="n";
    char line[]="cin";
    char copp[]="c";
    char unif[]="u";
    
    int nbbits=10;
    double alpha=1.4;
    int q;
    d=8;
    int m = 1;
    int sN,sX,h; 

    delta =0.75;

    int output = 0;

    PARSE_MAIN_ARGS {
      MATCH_MAIN_ARGID("-type",type);
      MATCH_MAIN_ARGID("-d",d);
      MATCH_MAIN_ARGID("-m",m);
      MATCH_MAIN_ARGID("-delta",delta);
      MATCH_MAIN_ARGID("-bits",nbbits);
      MATCH_MAIN_ARGID("-alpha",alpha);
      MATCH_MAIN_ARGID("-output",output);
      MATCH_MAIN_ARGID("-q",q);
      MATCH_MAIN_ARGID("-N",sN);
      MATCH_MAIN_ARGID("-X",sX);
      MATCH_MAIN_ARGID("-h",h);
    }

    // Knapsack 
    // --------
    if (strcmp(type,knapsack) ==0) {
     
      n=d+1;
      A.resize(d+1,d); 
      AT.resize(d,d+1);  
      AT.gen_intrel(nbbits);
      transpose(A,AT);
    } 

    // Block knapsack 
    // --------------
    else if (strcmp(type,block) ==0) {
      n=d;
      
      A.resize(d,d); 

      for (int i=0; i<m; i++)
	for (int j=0; j<d ; j++) {
	  A(i,j).randb(nbbits);

	}  
    
      Z_NR<mpz_t> one; 
      one = 1;
      for (int i=m; i<d; i++)
	A(i,i)=one;
    } 

    // Rectangular uniform
    // -------------------
    else if (strcmp(type,unif) ==0) {
      n=m;
      
      A.resize(n,d); 

      for (int i=0; i<n; i++)
	for (int j=0; j<d ; j++) {
	  A(i,j).randb(nbbits);

	}  
    } 

    // Ajtai 
    // -----
    else if (strcmp(type,ajtai) ==0) {
      n=d;
      A.resize(d,d); 
      AT.resize(d,d);  
      AT.gen_ajtai(alpha);
      transpose(A,AT);
    } 

    // Coppersmith PKC 14 BCFPNRZ 
    // -----
    else if (strcmp(type,copp) ==0) {

      int ddelta=3;
 
      n=ddelta*h;
      d=n;

      Z_NR<mpz_t> N,X,tz,one;
      one = 1;

      N.randb(sN);
      X.randb(sX);

      AT.resize(n,n);
      A.resize(n,n);

      int i,j,k;

      tz=one;
      for (i=h-1; i>=0; i--) {
	
	for (k=0; k<ddelta; k++) 
	  AT(i*ddelta+k,i*ddelta+k)=tz;
	
	tz.mul(tz,N);
      }

      tz=one;
      for (i=0; i<n; i++) {
	
	AT(i,i).mul(AT(i,i),tz);
	tz.mul(tz,X);
      }

      int s;
      
      for (i=0; i<n; i++) 
	for (j=i+1; j<n; j++) {
	  
	  s=size_in_bits(AT(i,i));
	  AT(j,i).randb(s);
	  
	}

      transpose(A,AT);

    } 

    // NTRU like 
    // ---------
    else if (strcmp(type,ntru) ==0) {
      d=2*d;
      n=d;
      A.resize(d,d); 
      AT.resize(d,d);  
      AT.gen_ntrulike(nbbits,q);
      transpose(A,AT);
    } 

    // -------------------------------
    else if (strcmp(type,line) ==0) {
      
      n=d;
      A.resize(d,d); 
      AT.resize(d,d); 

      cin >> AT ;
      
      transpose(A,AT);
    }

    else
      cout << "Warning: in function command_line_basis, no basis created" << endl; 

    
    if (output ==1) 
      print2maple(A,n,d);

   
  }



/* ***********************************************

          GENERATION   

   ********************************************** */

/* Using the defaulft mpfr precision */

/* A row real vector for one relation */
template<class RT> int gen3r2s(matrix<FP_NR<RT> >& B, int n, int r, int s) {

  B.resize(1,n); 

  FP_NR<RT> l2,l3; 

  FP_NR<RT> rr,ss;
 
  l2=2.0;
  l2.log(l2);
  ss=((double) s);   
  l2.div(l2,ss);
  l2.exponential(l2);

  l3=3.0;
  rr=((double) r);  
  l3.log(l3);
  l3.div(l3,rr);
  l3.exponential(l3);

  FP_NR<RT> alpha,beta;

  alpha.sub(l3,l2);
  alpha.log(alpha);

  beta=0.0;
  B(0,0)=1.0;
  for (int i=1; i<n; i++) {
    beta.add(beta,alpha);
    B(0,i).exponential(beta);
  }

  return 0;
}


/* Two row vectors for a simultaneaous relation */
template<class RT> int gen3r2s7t5u(matrix<FP_NR<RT> >& B, int n, int r, int s, int t, int u) {

  B.resize(2,n); 

  FP_NR<RT> l2,l3; 

  FP_NR<RT> rr,ss,tt,uu;

  // First vector 
  // ************
 
  l2=2.0;
  l2.log(l2);
  ss=((double) s);   
  l2.div(l2,ss);
  l2.exponential(l2);

  l3=3.0;
  rr=((double) r);  
  l3.log(l3);
  l3.div(l3,rr);
  l3.exponential(l3);

  FP_NR<RT> alpha,beta;

  alpha.sub(l3,l2);

  if (alpha  < 0) 
    alpha.abs(alpha);
  

  alpha.log(alpha);

  beta=0.0;
  B(0,0)=1.0;
  for (int i=1; i<n; i++) {
    beta.add(beta,alpha);
    B(0,i).exponential(beta);
    if ((i%2) > 0) B(0,i).neg(B(0,i));
  }

  // Second vector 
  // ************
 
  l2=5.0;
  l2.log(l2);
  uu=((double) u);   
  l2.div(l2,uu);
  l2.exponential(l2);

  l3=7.0;
  tt=((double) t);  
  l3.log(l3);
  l3.div(l3,tt);
  l3.exponential(l3);

  alpha.sub(l3,l2);

  if (alpha  < 0) 
    alpha.abs(alpha);
  

  alpha.log(alpha);

  beta=0.0;
  B(1,0)=1.0;
  for (int i=1; i<n; i++) {
    beta.add(beta,alpha);
    B(1,i).exponential(beta);
    if ((i%2) > 0) B(1,i).neg(B(1,i));
  }

  return 0;
}


} // end namespace hplll

#endif 
