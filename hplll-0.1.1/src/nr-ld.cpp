/*  long double FP_NR<long double> et ldpe_t 
    extra functions Z_NR and FP_NR 
    see nr.h and nr.cpp in the fplll library 
    hplllprint for parenthesis in print 
 
Copyright (C) 2011, 2012, 2013      Gilles Villard 

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

// Fichiers d'interface long double FP_NR<long double> et ldpe_t
// and extra functions, e.g. conversions 

#ifndef FPLLL_NR_LD
#define FPLLL_NR_LD

#include"defs.h"


FPLLL_BEGIN_NAMESPACE

// Rajouter si nul ou pas 
//***********************
// pas signe juste nul ou pas 

inline void get_round(FP_NR<dpe_t>& x, int&xs, long& lx, long& expo, const FP_NR<dpe_t> a, const FP_NR<dpe_t> b){

  DPE_DOUBLE qmant;
  DPE_EXP_T  qexp;

  qmant = DPE_MANT(a.getData()) / DPE_MANT(b.getData());
  qexp = DPE_EXP(a.getData()) - DPE_EXP(b.getData());

  // Normalization 
  DPE_EXP_T e;
  qmant = DPE_FREXP (qmant, &e);
  qexp += e;

  // Just zero, nothing assigned 
  if (qexp < 0) {
    xs = 0;
    return; 
  } 
  else {

    xs=1;

    if (qexp > DPE_BITSIZE) {

      if (sizeof(long) == 4) { /* 32-bit word: long has 31 bits */
	lx = (long) (qmant * 2147483648.0);
	expo = qexp - 31;
      }
      else if (sizeof(long) == 8) { /* 64-bit word: long has 63 bits */
	lx = (long) (qmant * 9223372036854775808.0);
	expo = qexp - 63;
      }	
      if (expo < 0 ) {
	lx = static_cast<long>(ldexp(static_cast<double>(lx), expo));
	expo=0;
      } 
      LDPE_MANT(x.getData())=qmant; 
      LDPE_EXP(x.getData())=qexp;
      

    } // End integer case   
    else {
      // sinon round via le double reformé ldexp puis remultiplié en long puis lexp
      // Un seul ldepx directement ? et cast double puis long cf get_si_exp fplll ? 

      qmant = round((ldexp(qmant, qexp)));
      lx = static_cast<long>(qmant);
      expo = 0; 

      LDPE_MANT(x.getData())=qmant; 
      LDPE_EXP(x.getData())=0;
      // Puis normaliser ? 
      qmant = DPE_FREXP (qmant, &e);
      LDPE_MANT(x.getData())=qmant; 
      LDPE_EXP(x.getData())=e;

     
    } 

  } // End non zero 

 

 


}


// Rajouté PSLQ passer de mpfr à ldpe
// ---------------------------------- 

inline void set_mpfr(FP_NR<dpe_t>& x, const FP_NR<mpfr_t> xf){

  long e;

  DPE_MANT(x.getData())=mpfr_get_d_2exp(&e, xf.getData(),GMP_RNDN); 
  DPE_EXP(x.getData())=e;
  dpe_normalize(x.getData());
}

#ifdef HPLLL_WITH_LONG_DOUBLE

inline void set_mpfr(FP_NR<ldpe_t>& x, const FP_NR<mpfr_t> xf){

  long e;

  LDPE_MANT(x.getData())=mpfr_get_ld_2exp(&e, xf.getData(),GMP_RNDN); 
  LDPE_EXP(x.getData())=e;
  ldpe_normalize(x.getData());
}
#endif 

// Strange ??? Some calls with mpfr in both matrices 
inline void set_mpfr(FP_NR<mpfr_t>& x, const FP_NR<mpfr_t> xf){

  x=xf;
}


inline void set_mpfr(FP_NR<double>& x, const FP_NR<mpfr_t> xf){

  x.getData() = xf.get_d();
}

#ifdef HPLLL_WITH_LONG_DOUBLE
inline void set_mpfr(FP_NR<long double>& x, const FP_NR<mpfr_t> xf){

  x.getData()= mpfr_get_ld(xf.getData(),GMP_RNDN);
}
#endif 


/* Rajouté pour ne pas toucher pour l'instant à l'interface nr.h */

template<class ZT, class FT> inline void set_f(Z_NR<ZT>& xz, const FP_NR<FT> x) {

  xz.set_f(x);

}

#ifdef HPLLL_WITH_LONG_DOUBLE

template<> inline void set_f(Z_NR<mpz_t>& xz, const FP_NR<long double> xx) {

  long double x;
  x=lround(xx.getData());

  int l;
  if (x==0) l=2;
  else if (x >0) 
    l=ilogbl(x)+2; // Longueur décimal bof
  else l=ilogbl(-x)+2; 

  char xc[l];    // À régler 

  sprintf(xc,"%.0Lf",x);

  mpz_set_str(xz.getData(),xc,10);

}
#endif 

#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void set_f(Z_NR<mpz_t>& xz, const FP_NR<ldpe_t> xx) {

  // pb dpe get_z en long double ?
  //ldpe_get_z(xz.getData(), const_cast<ldpe_t&>(x.GetData()));

  long double x;
  int l;

  if ( LDPE_EXP(xx.getData()) >= 64) { // We have an integer 

    x=ldexpl(LDPE_MANT(xx.getData()), 64); 
  
    if (x==0.0) l=2;
    else if (x > 0) 
      l=ilogbl(x)+2; // Longueur décimal bof
    else l=ilogbl(-x)+2; 

    char xc[l];    // À régler 

    sprintf(xc,"%.0Lf",x);
  
    mpz_set_str(xz.getData(),xc,10);

    mpz_mul_2exp (xz.getData(), xz.getData(), LDPE_EXP(xx.getData()) - 64);

  }
  else {

    x=ldexpl(LDPE_MANT(xx.getData()), LDPE_EXP(xx.getData()) ); 

    x=lround(x); // added April 10, 2013 - Needed ? Was required and missing in set_f from long double to mpz_t

    if (x==0.0) l=2;
    else if (x > 0) 
      l=ilogbl(x)+2; // Longueur décimal bof
    else l=ilogbl(-x)+2; 

    char xc[l];    // À régler 

    sprintf(xc,"%.0Lf",x);
  
    mpz_set_str(xz.getData(),xc,10);

  }  
}
#endif 

template<class ZT, class FT> inline void set_z(FP_NR<FT>& x, const Z_NR<ZT> xz) {

  x.set_z(xz);

}

#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void set_z(FP_NR<ldpe_t>& x, const Z_NR<mpz_t> xz){

  // cf pb set_z dpe en long double ? 
  //   ldpe_set_z(x.getData(), const_cast<mpz_t&>(xz.GetData()));

  mpfr_t xzf;

  mpfr_init2(xzf,80);  // À voir 

  long e;

  mpfr_set_z(xzf,xz.getData(),GMP_RNDN);

  LDPE_MANT(x.getData())=mpfr_get_ld_2exp(&e, xzf,GMP_RNDN); 
  LDPE_EXP(x.getData())=e;

  mpfr_clear(xzf);

}
#endif

#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void set_z(FP_NR<long double>& x, const Z_NR<mpz_t> xz){

  mpfr_t xzf;

  mpfr_init2(xzf,80);  // À voir 

  mpfr_set_z(xzf,xz.getData(),GMP_RNDN);

  x.getData()=mpfr_get_ld(xzf,GMP_RNDN); 

 mpfr_clear(xzf);

}
#endif 

// Bit size 
// --------

inline int size_in_bits(Z_NR<long> data) {
  long y = abs(data.GetData());
  int resul = 0;
  if (y == 0) resul=1;
  else {
    while (y > 1) {
      resul++;
      y >>= 1;  //y /= 2;
    }
  }
  return resul;
}

inline int size_in_bits(Z_NR<mpz_t> data) {
 
  int l;
  l= mpz_sizeinbase(data.getData(),2); 
  return l;
}





/**************** For print2maple **************/
// Put ^(-4) instead of ^-4 
//

template<class FT> inline void hplllprint(const FP_NR<FT> a);


template<> inline void hplllprint(const FP_NR<dpe_t> a) {


  DPE_DOUBLE d = DPE_MANT(a.getData());
  DPE_EXP_T e2 = DPE_EXP(a.getData());
  int e10 = 0;
  char sign = ' ';
  
  if (d == 0.0)  {
#ifdef DPE_USE_DOUBLE
    fprintf (stdout, "%1.*f", dpe_str_prec, d);
#else
    fprintf (stdout, "foo\n %1.*f", dpe_str_prec, d);
#endif
  }
  else { 

    if (d < 0)
      {
	d = -d;
	sign = '-';
      }

    if (e2 > 0)
      {
	while (e2 > 0)
	  {
	    e2 --;
	    d *= 2.0;
	    if (d >= 10.0)
	      {
		d /= 10.0;
		e10 ++;
	      }
	  }
      }
    else /* e2 <= 0 */
      {
	while (e2 < 0)
	  {
	    e2 ++;
	    d /= 2.0;
	    if (d < 1.0)
	      {
		d *= 10.0;
		e10 --;
	      }
	  }
      }

    
    fprintf (stdout, "%c%1.*f*10^(%d)", sign, dpe_str_prec, d, e10);
  }

  fflush(stdout);

}


#ifdef HPLLL_WITH_LONG_DOUBLE
template<> inline void hplllprint(const FP_NR<ldpe_t> a) {


  LDPE_DOUBLE d = LDPE_MANT(a.getData());
  LDPE_EXP_T e2 = LDPE_EXP(a.getData());
  int e10 = 0;
  char sign = ' ';
  
  if (d == 0.0)  {
#ifdef LDPE_USE_LONGDOUBLE
    fprintf (stdout, "%1.*Lf", ldpe_str_prec, d);
#else
    fprintf (stdout, "foo\n %1.*Lf", ldpe_str_prec, d);
#endif
  }
  else { 

    if (d < 0)
      {
	d = -d;
	sign = '-';
      }

    if (e2 > 0)
      {
	while (e2 > 0)
	  {
	    e2 --;
	    d *= 2.0;
	    if (d >= 10.0)
	      {
		d /= 10.0;
		e10 ++;
	      }
	  }
      }
    else /* e2 <= 0 */
      {
	while (e2 < 0)
	  {
	    e2 ++;
	    d /= 2.0;
	    if (d < 1.0)
	      {
		d *= 10.0;
		e10 --;
	      }
	  }
      }

    
    fprintf (stdout, "%c%1.*Lf*10^(%d)", sign, ldpe_str_prec, d, e10);
  }

  fflush(stdout);

}
#endif 

/*********
 * FP_NR *
 *********/
/*
template<class long double> inline void FP_NR<long double>::addmul(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t rnd)
{
  FP_NR<long double> product;
  product.mul(b, c, rnd);
  add(*this, product, rnd);
}

template<class F> inline void FP_NR<F>::submul(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd)
{
  FP_NR<F> product;
  product.mul(b, c, rnd);
  sub(*this, product, rnd);
}
*/

const int PREC_LONGDOUBLE=64;


/*************************
 * long double specialization *
 *************************/

// FROM FPLLL NOW some basics, not all  (April 2, 2013) 

#ifdef HPLLL_WITH_LONG_DOUBLE

template<> inline void FP_NR<long double>::set(const FP_NR<long double>& s)
{
  data=s.getData();
}


template<> inline void FP_NR<long double>::mul_2ui(const FP_NR<long double>& b, unsigned int c)
{
  data = ldexpl(b.GetData(), c);
}

template<> inline void FP_NR<long double>::div_2ui(const FP_NR<long double>& b, unsigned int c)
{
  data = ldexpl(b.GetData(), -c);
}


template<> inline int FP_NR<long double>::exp() const
{
  return ilogbl(data)+1;
}
template<> inline int FP_NR<long double>::zero_p() const
{
  return (data==0);
}

#endif 


/*************
 * LDPE spec. *
 *************/

#ifdef HPLLL_WITH_LONG_DOUBLE

#define ldpe_ncref const_cast<ldpe_t&>

template<> inline FP_NR<ldpe_t>::FP_NR()
{
  ldpe_init(data);
}
template<> inline FP_NR<ldpe_t>::FP_NR(const FP_NR<ldpe_t>& f)
{
  ldpe_init(data);
  ldpe_set(data, ldpe_ncref(f.data));
}
template<> inline FP_NR<ldpe_t>::~FP_NR()
{
  ldpe_clear(data);
}

template<> inline void FP_NR<ldpe_t>::print() const
{
  ldpe_out_str(stdout,10, ldpe_ncref(data));
  fflush(stdout);
}
template<> inline void FP_NR<ldpe_t>::printerr() const
{
  ldpe_out_str(stderr,10, ldpe_ncref(data));
  fflush(stderr);
}

template<> inline double FP_NR<ldpe_t>::get() const
{
  return ldpe_get_d(ldpe_ncref(data));
}
template<> inline double FP_NR<ldpe_t>::get_d(mp_rnd_t rnd) const
{
  return ldpe_get_d(ldpe_ncref(data));
}
template<> inline signed long int FP_NR<ldpe_t>::get_si() const
{
  return ldpe_get_si(ldpe_ncref(data));
}

// Added Lun 14 oct 2013 16:20:37 CEST / TO CHECK 
template<> 
inline long FP_NR<ldpe_t>::get_si_exp(long& expo) const {
  long result;
  expo = ldpe_get_si_exp(&result, ldpe_ncref(data));
  if (ldpe_zero_p(ldpe_ncref(data)))
    expo = 0;
  else if (expo < 0) {
    /* NOTE: conversion of result to double is exact even if
        sizeof(long) = 8 */
    result = static_cast<long>(ldexp(static_cast<long double>(result), expo));
    expo = 0;
  }
  return result;
}

template<> inline void FP_NR<ldpe_t>::set(const FP_NR<ldpe_t>& f)
{
  ldpe_set(data, ldpe_ncref(f.data));
}
template<> inline void FP_NR<ldpe_t>::set(double d)
{
  ldpe_set_d(data, d);
}
template<> inline void FP_NR<ldpe_t>::set(unsigned int s)
{
  ldpe_set_d(data, static_cast<double>(s));
}
template<> inline void FP_NR<ldpe_t>::operator=(const FP_NR<ldpe_t>& f)
{
  ldpe_set(data, ldpe_ncref(f.data));
}
template<> inline void FP_NR<ldpe_t>::operator=(double d)
{
  ldpe_set_d(data, d);
}

template<> inline void FP_NR<ldpe_t>::add(const FP_NR<ldpe_t>& a, const FP_NR<ldpe_t>& b, mp_rnd_t rnd)
{
  ldpe_add(data, ldpe_ncref(a.data), ldpe_ncref(b.data));
}
template<> inline void FP_NR<ldpe_t>::sub(const FP_NR<ldpe_t>& a, const FP_NR<ldpe_t>& b, mp_rnd_t rnd)
{
  ldpe_sub(data, ldpe_ncref(a.data), ldpe_ncref(b.data));
}
template<> inline void FP_NR<ldpe_t>::neg(const FP_NR<ldpe_t>& a)
{
  ldpe_neg(data, ldpe_ncref(a.data));
}
template<> inline void FP_NR<ldpe_t>::mul(const FP_NR<ldpe_t>& a, const FP_NR<ldpe_t>& b, mp_rnd_t rnd)
{
  ldpe_mul(data, ldpe_ncref(a.data), ldpe_ncref(b.data));
}
template<> inline void FP_NR<ldpe_t>::mul_2si(const FP_NR<ldpe_t>& a, long b)
{
  ldpe_mul_2si(data, ldpe_ncref(a.data),b);
}
template<> inline void FP_NR<ldpe_t>::mul_2ui(const FP_NR<ldpe_t>& a, unsigned int b)
{
  ldpe_mul_2exp(data, ldpe_ncref(a.data),b);
}
template<> inline void FP_NR<ldpe_t>::div(const FP_NR<ldpe_t>& a, const FP_NR<ldpe_t>& b, mp_rnd_t rnd)
{
  ldpe_div(data, ldpe_ncref(a.data), ldpe_ncref(b.data));
}
template<> inline void FP_NR<ldpe_t>::div_2ui(const FP_NR<ldpe_t>& a, unsigned int b)
{
  ldpe_div_2exp(data, ldpe_ncref(a.data),b);
}

template<> inline int FP_NR<ldpe_t>::cmp(const FP_NR<ldpe_t>& a) const
{
  return ldpe_cmp(ldpe_ncref(data), ldpe_ncref(a.data));
}
template<> inline int FP_NR<ldpe_t>::cmp(double a) const
{
  return ldpe_cmp_d(ldpe_ncref(data), a);
}
template<> inline int FP_NR<ldpe_t>::sgn() const
{
  return cmp(0.0);
}
template<> inline bool FP_NR<ldpe_t>::operator<=(const FP_NR<ldpe_t>& a) const
{
  return ldpe_cmp(ldpe_ncref(data), ldpe_ncref(a.data)) <= 0;
}
template<> inline bool FP_NR<ldpe_t>::operator<=(double a) const
{
  return ldpe_cmp_d(ldpe_ncref(data), a) <= 0;
}
template<> inline bool FP_NR<ldpe_t>::operator>=(const FP_NR<ldpe_t>& a) const
{
  return ldpe_cmp(ldpe_ncref(data), ldpe_ncref(a.data)) >= 0;
}
template<> inline bool FP_NR<ldpe_t>::operator>=(double a) const
{
  return ldpe_cmp_d(ldpe_ncref(data), a) >= 0;
}

template<> inline void FP_NR<ldpe_t>::abs(const FP_NR<ldpe_t>& a)
{
  ldpe_abs(data, ldpe_ncref(a.data));
}

template<> inline long FP_NR<ldpe_t>::exponent() const {
  return LDPE_EXP(data);
}

template<> inline void FP_NR<ldpe_t>::rnd(const FP_NR<ldpe_t>& a)
{
  ldpe_round(data, ldpe_ncref(a.data));
}
template<> inline int FP_NR<ldpe_t>::exp() const
{
  return LDPE_EXP(data);
}
template<> inline int FP_NR<ldpe_t>::zero_p() const
{
  return ldpe_zero_p(ldpe_ncref(data));
}
template<> inline void FP_NR<ldpe_t>::set_nan()
{
  //ldpe_set_d(data, NAN); // LDPE_UNLIKELY branch in ldpe_normalize
  LDPE_MANT(data) = NAN;
}
template<> inline int FP_NR<ldpe_t>::is_nan() const
{
  return (LDPE_MANT(data)!=LDPE_MANT(data));
}
template<> inline void FP_NR<ldpe_t>::sqrt(const FP_NR<ldpe_t>& a, mp_rnd_t rnd)
{
  ldpe_sqrt(data, ldpe_ncref(a.data));
}

template<>
inline void FP_NR<ldpe_t>::swap(FP_NR<ldpe_t>& a) {
  ldpe_swap(data, a.data);
}

template<>
inline unsigned int FP_NR<ldpe_t>::getprec() {
  return LDPE_BITSIZE;
}

template<>
inline unsigned int FP_NR<ldpe_t>::setprec(unsigned int prec) {
  // ignored
  return getprec();
}

#undef ldpe_ncref
#endif // LONG DOUBLE 

FPLLL_END_NAMESPACE

#endif

