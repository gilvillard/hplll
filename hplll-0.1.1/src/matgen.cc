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

  int sign;

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

  if (alpha  < 0) {
    sign=-1;
    alpha.abs(alpha);
  }

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

  if (alpha  < 0) {
    sign=-1;
    alpha.abs(alpha);
  }

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
