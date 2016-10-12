
/*********************************
 * Z=long 128 bits specialization
 *********************************/

#ifndef HPLLL_NR_Z_L128_H
#define FPLLL_NR_Z_L128_H

#include"defs.h"



FPLLL_BEGIN_NAMESPACE

// For gmp conversions 

// long on 64 bits 
// sign + 127 bits : 64 + 63 
inline void mpz_set_128int(Z_NR< __int128_t>& r, const Z_NR<mpz_t> a){

  mpz_t lmp,ump;
  mpz_init(lmp);
  mpz_init(ump);

  mpz_fdiv_q_2exp (ump, a.getData(), 64);
  mpz_fdiv_r_2exp (lmp, a.getData(), 64);

  r.getData()=mpz_get_si(ump);
  r.getData()<<=64;
  r.getData()+=mpz_get_ui(lmp);
}

// long on 64 bits 
// sign + 127 bits : 64 + 63 
inline void mpz_get_128int(Z_NR<mpz_t>& a, const Z_NR< __int128_t> r){

  mpz_t data;
  mpz_init(data);

  mpz_init(data);
  mpz_set_si(data,r.getData()>>64);
  mpz_mul_2exp(data,data,64);

  mpz_set(a.getData(),data);

  mpz_add_ui(a.getData(),a.getData(),(r.getData()<<64)>>64);

}


template<>
inline Z_NR<__int128_t>::Z_NR() {}

template<>
inline Z_NR<__int128_t>::Z_NR(const Z_NR<__int128_t>& z) : data(z.data) {}

template<>
inline Z_NR<__int128_t>::~Z_NR() {}


/** get data */
/* template<> */
/* inline double Z_NR<long>::get_d() const { */
/*   return static_cast<double>(data); */
/* } */

/* #ifdef FPLLL_WITH_LONG_DOUBLE */
/* template<> */
/* inline long double Z_NR<long>::get_ld() const { */
/*   return static_cast<long double>(data); */
/* } */
/* #endif */

/* template<> */
/* inline long Z_NR<long>::get_si() const { */
/*   return data; */
/* } */

/* template<> */
/* inline void Z_NR<long>::get_mpz(mpz_t r) const { */
/*   mpz_set_si (r, data); */
/* } */

/* inline long computeLongExponent(long data) { */
/*   unsigned long y = static_cast<unsigned long>(abs(data)); */
/*   long e; */
/*   for (e = 0; y; e++, y >>= 1) {} */
/*   return e; */
/* } */

/* template<> */
/* inline long Z_NR<long>::exponent() const { */
/*   int intExpo; */
/*   double fNorm = frexp(static_cast<double>(data), &intExpo); */
/*   if (data > MAX_LONG_FAST && fabs(fNorm) == 0.5) */
/*     return computeLongExponent(data); */
/*   else */
/*     return static_cast<long>(intExpo); */
/* } */

/* /\** set data *\/ */
/* template<> */
/* inline void Z_NR<long>::set_str(const char* s) { */
/*   data = atol(s); */
/* } */

/* /\** comparison *\/ */
template<>
inline int Z_NR<__int128_t>::cmp(const Z_NR<__int128_t>& m) const {
  if (data > m.data) return 1;
  if (data == m.data) return 0;
  return -1;
}

/* template<> */
/* inline int Z_NR<long>::sgn() const { */
/*   if (data > 0) return 1; */
/*   if (data == 0) return 0; */
/*   return -1; */
/* } */


/** operator */
template<>
inline void Z_NR<__int128_t>::operator=(const Z_NR<__int128_t>& a) {
  data = a.data;
}

template<>
inline void Z_NR<__int128_t>::operator=(const mpz_t& a) {
  data = static_cast<__int128_t>(mpz_get_si(a));
}

  

template<>
inline void Z_NR<__int128_t>::operator=(long a) {
  data = static_cast<__int128_t>(a);
}


template<> template<>
inline void Z_NR<__int128_t>::set_f(const FP_NR<double>& a) {
  //data = a.get_si();
  data = static_cast<__int128_t>(a.getData());
}


// Mettre ds nr-ld 
template<> template<>
inline void Z_NR<__int128_t>::set_f(const FP_NR<long double>& a) {
    data = static_cast<__int128_t>(a.getData());
}


template<> template<>
inline void Z_NR<__int128_t>::set_f(const FP_NR<dpe_t>& a) {
  data = a.get_si();
  //data = static_cast<__int128_t>(a.get_d());
}


//template<> template<>
//inline void Z_NR<__int128_t>::set_f(const FP_NR<mpfr_t>& a) {
//data = a.get_si();
//}
  
//template<> template<>
//inline void FP_NR<mpfr_t>::set_z(const Z_NR<__int128_t>& a, mp_rnd_t rnd) {
//mpfr_set_d(data, static_cast<double>(a.getData()), rnd);
//}

template<> template<>
inline void FP_NR<double>::set_z(const Z_NR<__int128_t>& a, mp_rnd_t rnd) {
  data=static_cast<double>(a.getData());
}

// Mettre ds nr-ld 
template<> template<>
inline void FP_NR<long double>::set_z(const Z_NR<__int128_t>& a, mp_rnd_t rnd) {
  data=static_cast<long double>(a.getData());
}




/* template<> */
/* inline void Z_NR<long>::operator=(const mpz_t& a) { */
/*   data = mpz_get_si(a); */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator<(const Z_NR<long>& a) const { */
/*   return data < a.data; */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator<(long a) const { */
/*   return data < a; */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator>(const Z_NR<long>& a) const { */
/*   return data > a.data; */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator>(long a) const { */
/*   return data > a; */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator<=(const Z_NR<long>& a) const { */
/*   return data <= a.data; */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator<=(long a) const { */
/*   return data <= a; */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator>=(const Z_NR<long>& a) const { */
/*   return data >= a.data; */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator>=(long a) const { */
/*   return data >= a; */
/* } */

template<>
inline bool Z_NR<__int128_t>::operator==(const Z_NR<__int128_t>& a) const {
  return data == a.data;
}

template<>
inline bool Z_NR<__int128_t>::operator==(long a) const {
  return data == a;
}

/* template<> */
/* inline bool Z_NR<long>::operator!=(const Z_NR<long>& a) const { */
/*   return data != a.data; */
/* } */

/* template<> */
/* inline bool Z_NR<long>::operator!=(long a) const { */
/*   return data != a; */
/* } */

/** arithmetic */
template<>
inline void Z_NR<__int128_t>::add(const Z_NR<__int128_t>& a, const Z_NR<__int128_t>& b) {
  data = a.data + b.data;
}

/* template<> */
/* inline void Z_NR<long>::add_ui(const Z_NR<long>& a, unsigned int b) { */
/*   data = a.data + b; */
/* } */

template<>
inline void Z_NR<__int128_t>::sub(const Z_NR<__int128_t>& a, const Z_NR<__int128_t>& b) {
  data = a.data - b.data;
}

/* template<> */
/* inline void Z_NR<long>::sub_ui(const Z_NR<long>& a, unsigned int b) { */
/*   data = a.data - b; */
/* } */

template<> 
inline void Z_NR<__int128_t>::neg(const Z_NR<__int128_t>& a) { 
  data = -a.data; 
} 

template<>
inline void Z_NR<__int128_t>::mul(const Z_NR<__int128_t>& a, const Z_NR<__int128_t>& b) {
  data = a.data * b.data;
}

template<>
inline void Z_NR<__int128_t>::mul_si(const Z_NR<__int128_t>& a, long b) {
  data = a.data * b;
}

/* template<> */
/* inline void Z_NR<long>::mul_ui(const Z_NR<long>& a, unsigned long b) { */
/*   data = a.data * b; */
/* } */

template<>
inline void Z_NR<__int128_t>::mul_2si(const Z_NR<__int128_t>& a, long b) {
  //NOTE: if b >= size_in_bits(a), the result is undefined
  if (b >= 0)
    data = a.data << b;
  else
    data = a.data >> -b;
}

/* template<> */
/* inline void Z_NR<long>::div_2si(const Z_NR<long>& a, long b) { */
/*   //NOTE: if b >= size_in_bits(a), the result is undefined */
/*   if (b >= 0) */
/*     data = a.data >> b; */
/*   else */
/*     data = a.data << -b; */
/* } */

template<>
inline void Z_NR<__int128_t>::addmul(const Z_NR<__int128_t>& a, const Z_NR<__int128_t>& b) {
  data += a.data * b.data;
}

/* template<> */
/* inline void Z_NR<long>::addmul_ui(const Z_NR<long>& a, unsigned long b) { */
/*   data += a.data * b; */
/* } */

template<>
inline void Z_NR<__int128_t>::addmul_si(const Z_NR<__int128_t>& a, long b) {
  data += a.data * b;
}

template<>
inline void Z_NR<__int128_t>::submul(const Z_NR<__int128_t>& a, const Z_NR<__int128_t>& b) {
  data -= a.data * b.data;
}

/* template<> */
/* inline void Z_NR<long>::submul_ui(const Z_NR<long>& a, unsigned long b) { */
/*   data -= a.data * b; */
/* } */

template<>
inline void Z_NR<__int128_t>::abs(const Z_NR<__int128_t>& a) {
  if (a.data >=0) 
    data = a.data;  //==> error Ven 13 mai 2016 18:34:13 CEST  with data = abs(a.data); __int128_t ???
  else data = -a.data;
}


/* template<> */
/* inline void Z_NR<long>::swap(Z_NR<long>& a) { */
/*   std::swap(data, a.data); */
/* } */

/** random numbers */
template<>
inline void Z_NR<__int128_t>::randb(int bits) {
  mpz_t temp;
  mpz_init(temp);
  mpz_urandomb(temp, RandGen::getGMPState(), bits);
  data = mpz_get_si(temp); 
  mpz_clear(temp);
}

/* template<> */
/* inline void Z_NR<long>::randb_si(int bits) { */
/*   randb (bits); */
/*   data = data * RandGenInt::getBit(); */
/* } */

/* template<> */
/* inline void Z_NR<long>::randm(const Z_NR<long>& max) { */
/*   mpz_t temp, lim; */
/*   mpz_init(temp); */
/*   mpz_init(lim); */
/*   mpz_set_si(lim, max.data); */
/*   mpz_urandomm(temp, RandGen::getGMPState(), lim); */
/*   data = mpz_get_si(temp); */
/*   mpz_clear(temp); */
/*   mpz_clear(lim); */
/* } */

/* template<> */
/* inline void Z_NR<long>::randm_si(const Z_NR<long>& max) { */
/*   randm (max); */
/*   data = data * RandGenInt::getBit(); */
/* } */


/* /\* FPLLL_V3_COMPAT *\/ */
/* #ifdef FPLLL_V3_COMPAT */

/* inline int sizeinbase2(long data) { */
/*   long y = abs(data); */
/*   int resul = 0; */
/*   if (y == 0) resul=1; */
/*   else { */
/*     while (y > 1) { */
/*       resul++; */
/*       y >>= 1;  //y /= 2; */
/*     } */
/*   } */
/*   return resul; */
/* } */

/* template<> */
/* inline void Z_NR<long>::print() const { */
/*   cout << data; */
/* } */

/* template<> */
/* inline void Z_NR<long>::printerr() const { */
/*   cerr << data; */
/* } */

/* template<> */
/* inline void Z_NR<long>::read() { */
/*   cin >> data; */
/* } */

/* template<> */
/* inline double Z_NR<long>::get_d_2exp(long* expo) const { */
/*   int intExpo = 0; */
/*   double x = frexp(static_cast<double>(data), &intExpo); */
/*   *expo = intExpo; */
/*   return x; */
/* } */

/* template<> */
/* inline void Z_NR<long>::set(/\*const*\/ long& s) { */
/*   data = s; */
/* } */

/* template<> */
/* inline void Z_NR<long>::set(unsigned long s) { */
/*   data = static_cast<long>(s); */
/* } */

/* #endif // #ifdef FPLLL_V3_COMPAT */


/* operators Z_NR<long> */
template<>
inline ostream& operator<<(ostream& os, const Z_NR<__int128_t>& x) {

__int128_t value= x.getData();

    std::ostream::sentry s( os );
    if ( s ) {
        __uint128_t tmp = value < 0 ? -value : value;
        char buffer[ 128 ];
        char* d = buffer +127;
        do
        {
            -- d;
            *d = "0123456789"[ tmp % 10 ];
            tmp /= 10;
        } while ( tmp != 0 );
        if ( value < 0 ) {
            -- d;
            *d = '-';
        }
        int len = buffer +127 - d;
        if ( os.rdbuf()->sputn( d, len ) != len ) {
            os.setstate( std::ios_base::badbit );
        }
    }
    return os;
}


/* template<> */
/* inline istream& operator>>(istream& is, Z_NR<__int128_t>& x) { */
/*   return is >> x.getData(); */
/* } */


FPLLL_END_NAMESPACE

#endif
