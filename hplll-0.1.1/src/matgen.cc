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
    char dec[]="dec";
    
    int nbbits=10;
    double alpha=1.4;
    int q;
    d=8;
    int m = 1;
    int sN,sX,h; 
    int shift =10;
    

    delta =0.99;

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
      MATCH_MAIN_ARGID("-shift",shift);
      
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

    // Upper part shifted 
    // ------------------
    else if (strcmp(type,dec) ==0) {

      n=d;      
      A.resize(d,d); 

      for (int i=1; i<d; i++)
	for (int j=0; j<d ; j++) {
	  A(i,j).randb(nbbits);
	}  

      for (int i=0; i<m; i++)
	for (int j=0; j<d; j++)
	  A(i,j).randb(nbbits+shift);

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
	  AT(j,i).randb(s+20);
	  
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

      for (int i=0; i<d/2; i++)
	for (int j=0; j<d; j++)
	  A(i,j)=AT(j,i+d/2);

      for (int i=0; i<d/2; i++)
	for (int j=0; j<d; j++)
	  A(i+d/2,j)=AT(j,i);
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

  
  // When testing prec and Seysen
  template<class ZT> int genalpha(ZZ_mat<ZT>& B, int n, double alpha) {

    int i,j;
    
    FP_NR<double> a0,ta;
    a0=alpha;
    
    B.resize(n,n);

    // Diagonal entries, decreasing powers of alpha
    // --------------------------------------------
    
    ta=a0;
    
    for (i=n-1; i>=0; i--) {  
      B(i,i).set_f(ta);
      ta.mul(ta,a0);
    }

    // Upper triangular entries, random less than the diagonal 
    // -------------------------------------------------------

    Z_NR<ZT> zrdm;
    
    FP_NR<mpfr_t> diagf,frdm, p49,p50;

    p49=1.0;
    p49.mul_2si(p49,49);
    
    p50=1.0;
    p50.mul_2si(p50,50);
   
    for (i=0; i<n; i++) {

      diagf.set_z(B(i,i));
		 
      for (j=i+1; j<n; j++) {

	zrdm.randb(50);
	
	frdm.set_z(zrdm);
	frdm.sub(frdm,p49);
	frdm.div(frdm,p50);	

	frdm.mul(frdm,diagf);
	
	B(i,j).set_f(frdm);

      }

    }


    for (i=0; i<n; i++) 
      for (j=0; j<i; j++) 
	B(i,j)=0;
	
    return 0;
  }
  
  
} // end namespace hplll

#endif 
