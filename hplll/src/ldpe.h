
/* Copyright (C) 2004, 2005, 2006, 2008 Patrick Pelissier, Paul Zimmermann,
  LORIA/INRIA Nancy - Grand-Est.

  GV 13-04-11  for long double only 

This file is part of the DPE Library.

The DPE Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The DPE Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the DPE Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

/* WARNING: Patched version */


using namespace std;

#include <stdlib.h> /* For abort */
#include <stdio.h>  /* For fprintf */
#include <math.h>   /* for round, floor, ceil */
#if defined (__sun) /* for round on Solaris 10 */ 
#include "tgmath.h" 
#endif 
#include <limits.h>

#ifndef __LDPE
#define __LDPE

#define LDPE_VERSION_MAJOR 1
#define LDPE_VERSION_MINOR 5

#if defined(__GNUC__) && (__GNUC__ >= 3)
# define LDPE_LIKELY(x) (__builtin_expect(!!(x),1))
# define LDPE_UNLIKELY(x) (__builtin_expect((x),0))
# define LDPE_UNUSED_ATTR  __attribute__((unused))
#else
# define LDPE_LIKELY(x) (x)
# define LDPE_UNLIKELY(x) (x)
# define LDPE_UNUSED_ATTR
#endif

#define LDPE_USE_LONGDOUBLE

#if defined(LDPE_USE_LONGDOUBLE)
# define LDPE_DOUBLE long double
# define LDPE_BITSIZE 64
# define LDPE_2_POW_BITSIZE 0x1P64
# define LDPE_LDEXP std::ldexp
# define LDPE_FREXP std::frexp
# define LDPE_ROUND roundl
# define LDPE_FLOOR std::floor
# define LDPE_CEIL std::ceil
# define LDPE_TRUNC truncl

#endif

#if defined(LDPE_USE_LONG)
# define LDPE_EXP_T  long    /* exponent type */
# define LDPE_EXPMIN LONG_MIN /* smallest possible exponent */
#elif defined(LDPE_USE_LONGLONG)
# define LDPE_EXP_T  long long
# define LDPE_EXPMIN LONG_LONG_MIN
#else
# define LDPE_EXP_T  int     /* exponent type */
# define LDPE_EXPMIN INT_MIN /* smallest possible exponent */
#endif

typedef union
{
 long double d;
 int i[2];
} ldpe_double_words;

typedef struct
{
 LDPE_DOUBLE d; /* significand */
 LDPE_EXP_T exp; /* exponent */
} ldpe_struct;

typedef ldpe_struct ldpe_t[1];

#define LDPE_MANT(x) ((x)->d)
#define LDPE_EXP(x)  ((x)->exp)
#define LDPE_SIGN(x) ((LDPE_MANT(x) < 0.0) ? -1 : (LDPE_MANT(x) > 0.0))

#define LDPE_INLINE static inline

/* initialize */
LDPE_INLINE void
ldpe_init (ldpe_t x LDPE_UNUSED_ATTR)
{
}

/* clear */
LDPE_INLINE void
ldpe_clear (ldpe_t x LDPE_UNUSED_ATTR)
{
}

/* set x to y */
LDPE_INLINE void
ldpe_set (ldpe_t x, ldpe_t y)
{
 LDPE_MANT(x) = LDPE_MANT(y);
 LDPE_EXP(x) = LDPE_EXP(y);
}

/* set x to -y */
LDPE_INLINE void
ldpe_neg (ldpe_t x, ldpe_t y)
{
 LDPE_MANT(x) = -LDPE_MANT(y);
 LDPE_EXP(x) = LDPE_EXP(y);
}

/* set x to |y| */
LDPE_INLINE void
ldpe_abs (ldpe_t x, ldpe_t y)
{
 LDPE_MANT(x) = (LDPE_MANT(y) >= 0) ? LDPE_MANT(y) : -LDPE_MANT(y);
 LDPE_EXP(x) = LDPE_EXP(y);
}

/* set mantissa in [1/2, 1), except for 0 which has minimum exponent */
/* FIXME: don't inline this function yet ? */
static void
ldpe_normalize (ldpe_t x)
{
 // GV Ven 21 fév 2014 10:13:24 CET
 //if (LDPE_UNLIKELY (LDPE_MANT(x) == 0.0 || finite (LDPE_MANT(x)) == 0))
 if (LDPE_UNLIKELY (LDPE_MANT(x) == 0.0 || isfinite (LDPE_MANT(x)) == 0))
   {
     if (LDPE_MANT(x) == 0.0)
       LDPE_EXP(x) = LDPE_EXPMIN;
     /* otherwise let the exponent of NaN, Inf unchanged */
   }
 else
   {
     LDPE_EXP_T e;

     long double m = LDPE_MANT(x);
     LDPE_MANT(x) = LDPE_FREXP (m, &e);
     LDPE_EXP(x) += e;

   }
}

LDPE_INLINE long double
ldpe_scale (long double d, int s)
{
 /* -LDPE_BITSIZE < s <= 0 and 1/2 <= d < 1 */

 return LDPE_LDEXP (d, s);

}

/* set x to y */
LDPE_INLINE void
ldpe_set_d (ldpe_t x, double y)
{
 LDPE_MANT(x) = (LDPE_DOUBLE) y;
 LDPE_EXP(x) = 0;
 ldpe_normalize (x);
}

/* set x to y */
LDPE_INLINE void
ldpe_set_ld (ldpe_t x, long double y)
{
 LDPE_MANT(x) = (LDPE_DOUBLE) y;
 LDPE_EXP(x) = 0;
 ldpe_normalize (x);
}

/* set x to y */
LDPE_INLINE void
ldpe_set_ui (ldpe_t x, unsigned long y)
{
 LDPE_MANT(x) = (LDPE_DOUBLE) y;
 LDPE_EXP(x) = 0;
 ldpe_normalize (x);
}

/* set x to y */
LDPE_INLINE void
ldpe_set_si (ldpe_t x, long y)
{
 LDPE_MANT(x) = (LDPE_DOUBLE) y;
 LDPE_EXP(x) = 0;
 ldpe_normalize (x);
}

LDPE_INLINE long
ldpe_get_si (ldpe_t x)
{
 LDPE_DOUBLE d = LDPE_LDEXP (LDPE_MANT (x), LDPE_EXP (x));
 return (long) d;
}

LDPE_INLINE unsigned long
ldpe_get_ui (ldpe_t x)
{
 LDPE_DOUBLE d = LDPE_LDEXP (LDPE_MANT (x), LDPE_EXP (x));
 return (d < 0.0) ? 0 : (unsigned long) d;
}

LDPE_INLINE double
ldpe_get_d (ldpe_t x)
{
 return LDPE_LDEXP (LDPE_MANT (x), LDPE_EXP (x));
}

LDPE_INLINE long double
ldpe_get_ld (ldpe_t x)
{
 return LDPE_LDEXP (LDPE_MANT (x), LDPE_EXP (x));
}

#ifdef __GMP_H__
/* set x to y */
// GV Attention pb en long double 
LDPE_INLINE void
ldpe_set_z (ldpe_t x, mpz_t y)
{
 long e;
 LDPE_MANT(x) = mpz_get_d_2exp (&e, y);
 LDPE_EXP(x) = (LDPE_EXP_T) e;
}

/* set x to y, rounded to nearest */
//GV Attention pb en long double 
LDPE_INLINE void
ldpe_get_z (mpz_t x, ldpe_t y)
{
 LDPE_EXP_T ey = LDPE_EXP(y);
 if (ey >= LDPE_BITSIZE) /* y is an integer */
   {
     LDPE_DOUBLE d = LDPE_MANT(y) * LDPE_2_POW_BITSIZE; /* d is an integer */
     mpz_set_d (x, d); /* should be exact */
     mpz_mul_2exp (x, x, (unsigned long) ey - LDPE_BITSIZE);
   }
 else /* LDPE_EXP(y) < LDPE_BITSIZE */
   {
     if (LDPE_UNLIKELY (ey < 0)) /* |y| < 1/2 */
       mpz_set_ui (x, 0);
     else
       {
         LDPE_DOUBLE d = LDPE_LDEXP(LDPE_MANT(y), ey);
         mpz_set_d (x, (long double) LDPE_ROUND(d));
       }
   }
}

/* return e and x such that y = x*2^e */
LDPE_INLINE mp_exp_t
ldpe_get_z_exp (mpz_t x, ldpe_t y)
{
 mpz_set_d (x, LDPE_MANT (y) * LDPE_2_POW_BITSIZE);
 return LDPE_EXP(y) - LDPE_BITSIZE;
}
#endif

/* x <- y + z, assuming y and z are normalized, returns x normalized */
LDPE_INLINE void
ldpe_add (ldpe_t x, ldpe_t y, ldpe_t z)
{
 if (LDPE_UNLIKELY (LDPE_EXP(y) > LDPE_EXP(z) + LDPE_BITSIZE))
   /* |z| < 1/2*ulp(y), thus o(y+z) = y */
   ldpe_set (x, y);
 else if (LDPE_UNLIKELY (LDPE_EXP(z) > LDPE_EXP(y) + LDPE_BITSIZE))
   ldpe_set (x, z);
 else
   {
     LDPE_EXP_T d = LDPE_EXP(y) - LDPE_EXP(z); /* |d| <= LDPE_BITSIZE */

     if (d >= 0)
       {
         LDPE_MANT(x) = LDPE_MANT(y) + ldpe_scale (LDPE_MANT(z), -d);
         LDPE_EXP(x) = LDPE_EXP(y);
       }
     else
       {
         LDPE_MANT(x) = LDPE_MANT(z) + ldpe_scale (LDPE_MANT(y), d);
         LDPE_EXP(x) = LDPE_EXP(z);
       }
     ldpe_normalize (x);
   }
}

/* x <- y - z, assuming y and z are normalized, returns x normalized */
LDPE_INLINE void
ldpe_sub (ldpe_t x, ldpe_t y, ldpe_t z)
{
 if (LDPE_UNLIKELY (LDPE_EXP(y) > LDPE_EXP(z) + LDPE_BITSIZE))
   /* |z| < 1/2*ulp(y), thus o(y-z) = y */
   ldpe_set (x, y);
 else if (LDPE_UNLIKELY (LDPE_EXP(z) > LDPE_EXP(y) + LDPE_BITSIZE))
   ldpe_neg (x, z);
 else
   {
     LDPE_EXP_T d = LDPE_EXP(y) - LDPE_EXP(z); /* |d| <= LDPE_BITSIZE */

     if (d >= 0)
       {
         LDPE_MANT(x) = LDPE_MANT(y) - ldpe_scale (LDPE_MANT(z), -d);
         LDPE_EXP(x) = LDPE_EXP(y);
       }
     else
       {
         LDPE_MANT(x) = ldpe_scale (LDPE_MANT(y), d) - LDPE_MANT(z);
         LDPE_EXP(x) = LDPE_EXP(z);
       }
     ldpe_normalize (x);
   }
}

/* x <- y * z, assuming y and z are normalized, returns x normalized */
LDPE_INLINE void
ldpe_mul (ldpe_t x, ldpe_t y, ldpe_t z)
{
 LDPE_MANT(x) = LDPE_MANT(y) * LDPE_MANT(z);
 LDPE_EXP(x) = LDPE_EXP(y) + LDPE_EXP(z);
 ldpe_normalize (x);
}

/* x <- sqrt(y), assuming y is normalized, returns x normalized */
LDPE_INLINE void
ldpe_sqrt (ldpe_t x, ldpe_t y)
{
 LDPE_EXP_T ey = LDPE_EXP(y);
 if (ey % 2)
   {
     /* since 1/2 <= my < 1, 1/4 <= my/2 < 1 */
     LDPE_MANT(x) = std::sqrt(0.5 * LDPE_MANT(y));
     LDPE_EXP(x) = (ey + 1) / 2;
   }
 else
   {
     LDPE_MANT(x) = std::sqrt(LDPE_MANT(y));
     LDPE_EXP(x) = ey / 2;
   }
}

/* x <- y / z, assuming y and z are normalized, returns x normalized.
  Assumes z is not zero. */
LDPE_INLINE void
ldpe_div (ldpe_t x, ldpe_t y, ldpe_t z)
{
 LDPE_MANT(x) = LDPE_MANT(y) / LDPE_MANT(z);
 LDPE_EXP(x) = LDPE_EXP(y) - LDPE_EXP(z);
 ldpe_normalize (x);
}

/* x <- y * z, assuming y normalized, returns x normalized */
LDPE_INLINE void
ldpe_mul_ui (ldpe_t x, ldpe_t y, unsigned long z)
{
 LDPE_MANT(x) = LDPE_MANT(y) * (LDPE_DOUBLE) z;
 LDPE_EXP(x) = LDPE_EXP(y);
 ldpe_normalize (x);
}

/* x <- y / z, assuming y normalized, z non-zero, returns x normalized */
LDPE_INLINE void
ldpe_div_ui (ldpe_t x, ldpe_t y, unsigned long z)
{
 LDPE_MANT(x) = LDPE_MANT(y) / (LDPE_DOUBLE) z;
 LDPE_EXP(x) = LDPE_EXP(y);
 ldpe_normalize (x);
}

/* x <- y * 2^e */
LDPE_INLINE void
ldpe_mul_2si (ldpe_t x, const ldpe_t y, long e)
{
  LDPE_MANT(x) = LDPE_MANT(y);
  LDPE_EXP(x) = LDPE_EXP(y) + (LDPE_EXP_T) e;
}

/* x <- y * 2^e */
LDPE_INLINE void
ldpe_mul_2exp (ldpe_t x, ldpe_t y, unsigned long e)
{
 LDPE_MANT(x) = LDPE_MANT(y);
 LDPE_EXP(x) = LDPE_EXP(y) + (LDPE_EXP_T) e;
}

/* x <- y / 2^e */
LDPE_INLINE void
ldpe_div_2exp (ldpe_t x, ldpe_t y, unsigned long e)
{
 LDPE_MANT(x) = LDPE_MANT(y);
 LDPE_EXP(x) = LDPE_EXP(y) - (LDPE_EXP_T) e;
}

/* return e and x such that y = x*2^e (equality is not guaranteed if the 'long'
  type has fewer bits than the significand in ldpe_t) */
LDPE_INLINE mp_exp_t
ldpe_get_si_exp (long *x, ldpe_t y)
{
 if (sizeof(long) == 4) /* 32-bit word: long has 31 bits */
   {
     *x = (long) (LDPE_MANT(y) * 2147483648.0);
     return LDPE_EXP(y) - 31;
   }
 else if (sizeof(long) == 8) /* 64-bit word: long has 63 bits */
   {
     *x = (long) (LDPE_MANT (y) * 9223372036854775808.0);
     return LDPE_EXP(y) - 63;
   }
 else
   {
     fprintf (stderr, "Error, neither 32-bit nor 64-bit word\n");
     exit (1);
   }
}

static LDPE_UNUSED_ATTR int ldpe_str_prec = 16;
static int ldpe_out_str (FILE *s, int base, ldpe_t x) LDPE_UNUSED_ATTR;

static int
ldpe_out_str (FILE *s, int base, ldpe_t x)
{
 LDPE_DOUBLE d = LDPE_MANT(x);
 LDPE_EXP_T e2 = LDPE_EXP(x);
 int e10 = 0;
 char sign = ' ';
 if (LDPE_UNLIKELY (base != 10))
   {
     fprintf (stderr, "Error in ldpe_out_str, only base 10 allowed\n");
     exit (1);
   }
 if (d == 0.0)
#ifdef LDPE_USE_LONGDOUBLE
   return fprintf (s, "%1.*Lf", dpe_str_prec, d);
#else
   return fprintf (s, "foo\n %1.*Lf", dpe_str_prec, d);
#endif

   //return fprintf (s, "foo\n %1.*Lf", ldpe_str_prec, d);
   //return fprintf (s, "%1.*Lf", ldpe_str_prec, d);

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

 return fprintf (s, "%c%1.*Lf*10^%d", sign, ldpe_str_prec, d, e10);

}

static size_t ldpe_inp_str (ldpe_t x, FILE *s, int base) LDPE_UNUSED_ATTR;

static size_t
ldpe_inp_str (ldpe_t x, FILE *s, int base)
{
 size_t res;
 LDPE_DOUBLE d;
 if (LDPE_UNLIKELY (base != 10))
   {
     fprintf (stderr, "Error in ldpe_out_str, only base 10 allowed\n");
     exit (1);
   }

 res = fscanf (s, "%Lf", &d);

 ldpe_set_d (x, d);
 return res;
}

LDPE_INLINE void
ldpe_dump (ldpe_t x)
{
 ldpe_out_str (stdout, 10, x);
 putchar ('\n');
}

LDPE_INLINE int
ldpe_zero_p (ldpe_t x)
{
 return LDPE_MANT (x) == 0;
}

/* return a positive value if x > y
         a negative value if x < y
         and 0 otherwise (x=y). */
LDPE_INLINE int
ldpe_cmp (ldpe_t x, ldpe_t y)
{
 int sx = LDPE_SIGN(x);
 int d = sx - LDPE_SIGN(y);

 if (d != 0)
   return d;
 else if (LDPE_EXP(x) > LDPE_EXP(y))
   return (sx > 0) ? 1 : -1;
 else if (LDPE_EXP(y) > LDPE_EXP(x))
   return (sx > 0) ? -1 : 1;
 else /* LDPE_EXP(x) = LDPE_EXP(y) */
   return (LDPE_MANT(x) < LDPE_MANT(y)) ? -1 : (LDPE_MANT(x) > LDPE_MANT(y));
}

LDPE_INLINE int
ldpe_cmp_d (ldpe_t x, long double d)
{
 ldpe_t y;
 ldpe_set_d (y, d);
 return ldpe_cmp (x, y);
}

LDPE_INLINE int
ldpe_cmp_ui (ldpe_t x, unsigned long d)
{
 ldpe_t y;
 ldpe_set_ui (y, d);
 return ldpe_cmp (x, y);
}

LDPE_INLINE int
ldpe_cmp_si (ldpe_t x, long d)
{
 ldpe_t y;
 ldpe_set_si (y, d);
 return ldpe_cmp (x, y);
}

/* set x to integer nearest to y */
LDPE_INLINE void
ldpe_round (ldpe_t x, ldpe_t y)
{
 if (LDPE_EXP(y) < 0) /* |y| < 1/2 */
   ldpe_set_ui (x, 0);
 else if (LDPE_EXP(y) >= LDPE_BITSIZE) /* y is an integer */
   ldpe_set (x, y);
 else
   {
     LDPE_DOUBLE d;
     d = LDPE_LDEXP(LDPE_MANT(y), LDPE_EXP(y));
     ldpe_set_d (x, LDPE_ROUND(d));
   }
}

/* set x to the fractional part of y, defined as y - trunc(y), thus the
  fractional part has absolute value in [0, 1), and same sign as y */
LDPE_INLINE void
ldpe_frac (ldpe_t x, ldpe_t y)
{
 /* If |y| is smaller than 1, keep it */
 if (LDPE_EXP(y) <= 0)
   ldpe_set (x, y);
 else if (LDPE_EXP(y) >= LDPE_BITSIZE) /* y is an integer */
   ldpe_set_ui (x, 0);
 else
   {
     LDPE_DOUBLE d;
     d = LDPE_LDEXP(LDPE_MANT(y), LDPE_EXP(y));
     ldpe_set_d (x, d - LDPE_TRUNC(d));
   }
}

/* set x to largest integer <= y */
LDPE_INLINE void
ldpe_floor (ldpe_t x, ldpe_t y)
{
 if (LDPE_EXP(y) <= 0) /* |y| < 1 */
   {
     if (LDPE_SIGN(y) >= 0) /* 0 <= y < 1 */
       ldpe_set_ui (x, 0);
     else /* -1 < y < 0 */
       ldpe_set_si (x, -1);
   }
 else if (LDPE_EXP(y) >= LDPE_BITSIZE) /* y is an integer */
   ldpe_set (x, y);
 else
   {
     LDPE_DOUBLE d;
     d = LDPE_LDEXP(LDPE_MANT(y), LDPE_EXP(y));
     ldpe_set_d (x, LDPE_FLOOR(d));
   }
}

/* set x to smallest integer >= y */
LDPE_INLINE void
ldpe_ceil (ldpe_t x, ldpe_t y)
{
 if (LDPE_EXP(y) <= 0) /* |y| < 1 */
   {
     if (LDPE_SIGN(y) > 0) /* 0 < y < 1 */
       ldpe_set_ui (x, 1);
     else /* -1 < y <= 0 */
       ldpe_set_si (x, 0);
   }
 else if (LDPE_EXP(y) >= LDPE_BITSIZE) /* y is an integer */
   ldpe_set (x, y);
 else
   {
     LDPE_DOUBLE d;
     d = LDPE_LDEXP(LDPE_MANT(y), LDPE_EXP(y));
     ldpe_set_d (x, LDPE_CEIL(d));
   }
}

LDPE_INLINE void
ldpe_swap (ldpe_t x, ldpe_t y)
{
 LDPE_EXP_T i = LDPE_EXP (x);
 LDPE_DOUBLE d = LDPE_MANT (x);
 LDPE_EXP (x) = LDPE_EXP (y);
 LDPE_MANT (x) = LDPE_MANT (y);
 LDPE_EXP (y) = i;
 LDPE_MANT (y) = d;
}

#endif /* __LDPE */


