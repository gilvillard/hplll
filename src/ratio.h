/* Basis quality 

Created Mer  3 fév 2016 15:47:55 CET
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


#ifndef HPLLL_RATIO_H
#define HPLLL_RATIO_H

namespace hplll {


// ******************************************************************************
//
//  Statitics concerning the quality of the reduction
//  The matrix is assumed to be reduced   n rows >= d columns
//
//   log_2 of the Frobenius norm cond 
//
// ******************************************************************************


template<class ZT> int ratio(ZZ_mat<ZT> B, double& lfcond,  double& av_ratio,  double& max_ratio, double& c) {

  int i,j,k;


  //int n=B.getRows();
  int d=B.get_cols();
  
  mpfr_set_default_prec(2*d);
  
  Lattice<ZT, mpfr_t, matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > L(B);

  L.householder();

  matrix<FP_NR<mpfr_t> > R;

  R= L.getR();

  // ** Frobenius norm cond  ** 
  // **************************

  matrix<FP_NR<mpfr_t> > aR;
  matrix<FP_NR<mpfr_t> > iR;

  // Absolute value
  
  aR.resize(d,d);
      
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      aR(i,j).abs(R(i,j));

  // Inversion 
  
  iR.resize(d,d);
  
  for (i=0; i<d; i++) iR(i,i)=1.0; 
  
  // Higham p263 Method 1
  
  FP_NR<mpfr_t> t,one;
  one=1.0;

  for (j=0; j<d; j++) {
    iR(j,j).div(one,aR(j,j));
    for (k=j+1; k<d; k++) {
      t.mul(iR(j,j),R(j,k));
      iR(j,k).neg(t);
    }
    
    for (k=j+1; k<d; k++) {
      for (i=j+1; i<k; i++) {
	t.mul(iR(j,i),R(i,k)); 
	iR(j,k).sub(iR(j,k),t); 
      }
      iR(j,k).div(iR(j,k),R(k,k)); 
    }
  }

  // Abs iR --> iR 
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      iR(i,j).abs(iR(i,j));
  
  // aR * iR 
  
  matrix<FP_NR<mpfr_t> > prod;
  prod.resize(d,d);
  
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) {
      prod(i,j).mul(aR(i,0),iR(0,j));
      for (k=1; k<d ; k++) 
	prod(i,j).addmul(aR(i,k),iR(k,j));
    }


  // Frobenius norm 

  FP_NR<mpfr_t>  cc;
  cc=0.0; 
  
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      cc.addmul(prod(i,j),prod(i,j));
  
  cc.sqrt(cc); 

  mpfr_log2(cc.get_data(),cc.get_data(),GMP_RNDN);

  lfcond = cc.get_d();

  // ** Diagonal ratio ** 
  // ********************

  FP_NR<mpfr_t>  maxq,s,tt;

  s=0.0;
  maxq = 0.0;

  for (i=1; i<d; i++) {

    tt.div(aR(i-1,i-1),aR(i,i));

    s.add(s,tt);

    if (tt.cmp(maxq) > 0) maxq=tt;

  }

  // Average diagonal ratio
  
  tt= ((double) d);
  s.div(s,tt);

  av_ratio = s.get_d();
  max_ratio = maxq.get_d();


  // ** First vector length ** 
  // *************************


  // First vector 
  mpfr_log2(t.get_data(),(aR(0,0)).get_data(),GMP_RNDN);
  
  // Volume
  
  FP_NR<mpfr_t>  v;

  v=aR(0,0);
  for (i=1; i<d; i++) v.mul(v,aR(i,i));
  mpfr_log2(v.get_data(),v.get_data(),GMP_RNDN);

  tt= ((double) d);
  v.div(v,tt);

  t.sub(t,v);
  t.div(t,tt);

  // Or keep t?
  c = pow(2.0,t.get_d());
  
  return 0;
  
}

} // end namespace hplll


#endif


