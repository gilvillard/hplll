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


using namespace std;

#include"defs.h"

FPLLL_BEGIN_NAMESPACE

// Rajouté PSLQ passer de mpfr à ldpe
// ---------------------------------- 

inline void set_mpfr(FP_NR<dpe_t>& x, const FP_NR<mpfr_t> xf){

  long e;

  DPE_MANT(x.getData())=mpfr_get_d_2exp(&e, xf.getData(),GMP_RNDN); 
  DPE_EXP(x.getData())=e;
  dpe_normalize(x.getData());
}

inline void set_mpfr(FP_NR<ldpe_t>& x, const FP_NR<mpfr_t> xf){

  long e;

  LDPE_MANT(x.getData())=mpfr_get_ld_2exp(&e, xf.getData(),GMP_RNDN); 
  LDPE_EXP(x.getData())=e;
  ldpe_normalize(x.getData());
}

// Stange ??? Some calls with mpfr in both matrices 
inline void set_mpfr(FP_NR<mpfr_t>& x, const FP_NR<mpfr_t> xf){

  x=xf;
}


inline void set_mpfr(FP_NR<double>& x, const FP_NR<mpfr_t> xf){

  x.getData() = xf.get_d();
}

inline void set_mpfr(FP_NR<long double>& x, const FP_NR<mpfr_t> xf){

  x.getData()= mpfr_get_ld(xf.getData(),GMP_RNDN);
}



/* Rajouté pour ne pas toucher pour l'instant à l'interface nr.h */

template<class ZT, class FT> inline void set_f(Z_NR<ZT>& xz, const FP_NR<FT> x) {

  xz.set_f(x);

}

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


template<class ZT, class FT> inline void set_z(FP_NR<FT>& x, const Z_NR<ZT> xz) {

  x.set_z(xz);

}


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

template<> inline void set_z(FP_NR<long double>& x, const Z_NR<mpz_t> xz){

  mpfr_t xzf;

  mpfr_init2(xzf,80);  // À voir 

  mpfr_set_z(xzf,xz.getData(),GMP_RNDN);

  x.getData()=mpfr_get_ld(xzf,GMP_RNDN); 

 mpfr_clear(xzf);

}


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

// FROM FPLLL NOW (April 2, 2013) 
#ifndef FPLLL_WITH_LONG_DOUBLE

template<> inline FP_NR<long double>::FP_NR()
{
}

template<> inline FP_NR<long double>::FP_NR(const FP_NR<long double>& f)
{
  data = f.data;
}


template<> inline FP_NR<long double>::~FP_NR()
{
}


template<> inline void FP_NR<long double>::print() const
{
  cout << data;
}
template<> inline void FP_NR<long double>::printerr() const
{
  cerr << data;
}


template<> inline double FP_NR<long double>::get_d(mp_rnd_t rnd) const
{
  return static_cast<long double>(data);
  }
template<> inline signed long int FP_NR<long double>::get_si() const
{
  return static_cast<signed long int>(data);
}

template<> inline void FP_NR<long double>::operator=(const FP_NR<long double>& f)
{
  data = f.data;
}
template<> inline void FP_NR<long double>::operator=(double d)
{
  data = static_cast<long double>(d);
}


template<> inline void FP_NR<long double>::add(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t rnd)
{
  data= b.getData()+c.getData();
}
template<> inline void FP_NR<long double>::sub(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t rnd)
{
  data= b.getData()-c.getData();
}
template<> inline void FP_NR<long double>::neg(const FP_NR<long double>& b)
{
  data=-b.GetData();
}
template<> inline void FP_NR<long double>::mul(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t rnd)
{
  data= b.GetData()*c.GetData();
}

template<> inline void FP_NR<long double>::div(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t rnd)
{
  data= b.GetData()/c.GetData();
}

template<> inline int  FP_NR<long double>::cmp(const FP_NR<long double>& b) const {
  if (data < b.data)
    return -1;
  else if (data > b.data)
    return 1;
  else
    return 0;
}


template<> inline int  FP_NR<long double>::cmp(double b) const
{
  if (data<b)
    return -1;
  else if (data > b) 
    return 1;
  else 
    return 0;
}
template<> inline int FP_NR<long double>::sgn() const
{
  if (data>0)
    return 1;
  else if (data==0)
    return 0;
  else return -1;
  }
template<> inline bool FP_NR<long double>::operator<=(const FP_NR<long double>& a) const
{
  return data <= a.data;
}
template<> inline bool FP_NR<long double>::operator<=(double a) const
{
  return data <= a;
  }
template<> inline bool FP_NR<long double>::operator>=(const FP_NR<long double>& a) const
{
  return data >= a.data;
}

template<> inline bool FP_NR<long double>::operator>=(double a) const
{
  return data >= a;
  }


template<> inline void FP_NR<long double>::abs(const FP_NR<long double>& b)
{
  data=b.GetData();
  if (data<0) data=-data;
}
template<> inline void FP_NR<long double>::rnd(const FP_NR<long double>& b)
{
  data=rintl(b.GetData());
}

template<> inline void FP_NR<long double>::set_nan()
{
  data=NAN; // 0.0/0.0;
}
template<> inline int FP_NR<long double>::is_nan() const
{
  return (data!=data);
}
template<> inline void FP_NR<long double>::sqrt(const FP_NR<long double>& s, mp_rnd_t rnd)
{
  data = std::sqrt(s.GetData());  // sqrtl
}

template<>
inline void FP_NR<long double>::pow_si(const FP_NR<long double>& a, long int b, mp_rnd_t rnd) {
  data = std::pow(a.data, static_cast<long double>(b));
}

template<> inline void FP_NR<long double>::exponential(const FP_NR<long double>& a, mp_rnd_t rnd)
{
  data = std::exp(a.data);
}
template<> inline void FP_NR<long double>::log(const FP_NR<long double>& a, mp_rnd_t rnd)
{
  data = std::log(a.data);
}

template<>
inline void FP_NR<long double>::swap(FP_NR<long double>& a) {
  std::swap(data, a.data);
}


template<>
inline unsigned int FP_NR<long double>::getprec() {
  return PREC_LONGDOUBLE;
}

template<>
inline unsigned int FP_NR<long double>::setprec(unsigned int prec) {
  // ignored
  return getprec();
}

#endif  // From FPLLL 


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




/*************
 * LDPE spec. *
 *************/

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

template<>
inline long FP_NR<ldpe_t>::exponent() const {
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

FPLLL_END_NAMESPACE

#endif

