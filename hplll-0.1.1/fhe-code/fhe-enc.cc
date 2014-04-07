// Copyright (C) IBM, All Rights Reserved

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

NTL_CLIENT
#include "randomStuff.h"
#include "fhe.h"

static void evalRandPoly(vec_ZZ& vals, ZZ& r_to_m,
			 long n, long d, double p, const ZZ& r, const ZZ& M);
static void basicRandPoly(vec_ZZ& vals, ZZ& r_to_m,
			  long n, long d, double p, const ZZ& r, const ZZ& M);

/** Batch encryption of num bits. In principle, the encryption
 * algorithm, with pubkey (det,root) and bit b, is:
 *   1. Choose a random polynomial u() with one-bit coefficients
 *   2. Set the polynomial a() as a() = 2*u() + b
 *   3. Output the integer a(root) mod det
 ********************************************************************/
bool FHEkeys::encrypt(vec_ZZ& c, unsigned int b[], int num) const
{
  int i;
  ZZ tmp;       // used to hold r^m in the recursive calls
  unsigned long n = 1UL<<(prms.logn); // the dimenssion
  double p = ((double)prms.noise)/n;  // # of expected nonzero coefficients
  if (p>0.5) p = 0.5;

  // Evaluate all the polynomials together at root mode det
  evalRandPoly(c, tmp, num, n, p, root, det);

  // Set c[i] = 2*c[i] + b[i]
  for (i=0; i<num; i++) {
    c[i] <<= 1;
    c[i] += b[i];
    if (c[i]>=det) c[i] -= det;
  }
  return true;
}

/********************************************************************/
/******                   Utility functions                    ******/
/********************************************************************/

// Choose n random degree-(m-1) 0/1 polynomials, and evaluate them at r mod M.
// Returns the results in the vector vals, and return r^m in r_to_m. Upto
// a threshold it just calls the "trivial procedure" basicRandPoly, and
// above this threshold it recursively calls evalRandPoly(...,2n,m/2,...)

inline bool aboveThreshold(long n, long m) { return (n+1+(m&1)<m/2); }
static void evalRandPoly(vec_ZZ& vals, ZZ& r_to_m,
			 long n, long m, double p, const ZZ& r, const ZZ& M)
{
  if (aboveThreshold(n,m)) { // make a recursive call
    evalRandPoly(vals, r_to_m, 2*n, m/2, p, r, M); // r_to_m = r^{m/2}
    ZZ tmp;
    int s;
    for (int i=0; i<n; i++) { // val[i] += r^{m/2} * val[i+n] mod M
      // If m is odd then add another random 0/1 coefficient
      if ((m&1) && (s=randomBit(p)))
	if (s==1) AddMod(vals[i+n],vals[i+n],r_to_m,M);
	else      SubMod(vals[i+n],vals[i+n],r_to_m,M);
      MulMod(tmp, vals[i+n], r_to_m, M);// multiply "top half" by r^{d/2}
      vals[i] += tmp;
    }
    // free unneeded space
    for (int i=0; i<n; i++) vals[i+n].kill();
    vals.SetLength(n);

    SqrMod(r_to_m, r_to_m, M);          // compute r^m for the next level
    if (m&1) MulMod(r_to_m,r_to_m,r,M); // if m is odd, multiply by r again
  }
  else basicRandPoly(vals, r_to_m, n, m, p, r, M);
}


// Choose n random degree-m 0/1 polynomials, and evaluate them at r mod M.
// Returns the results in the vector vals, and return r^m in r_to_m. This
// procedure uses the "trivial" algorithm that just copmutes all the powers
// r,r^2,...,r^m and adds the corresponding power for each polynomial.
// This implementation uses modular squaring (rather than multiplication)
// to compute the even powers r^{2j}
static void basicRandPoly(vec_ZZ& vals, ZZ& r_to_m,
			  long n, long m, double p, const ZZ& r, const ZZ& M)
{
  int i,j,k,s;
  vals.SetLength(n);
  if (m<=0) {
    r_to_m = to_ZZ(1);
    return;
  }

  for (i=0; i<n; i++) vals[i] = randomBit(p); // the free term (0/1)
  if (m==1) {
    r_to_m = r;
    return;
  }
    
  ZZ rSqr; SqrMod(rSqr,r,M); // holds the value r^2 mod M

  // Handle the powers 1,2,4,... separately (saves maybe 1-2 mults)
  for (i=0; i<n; i++) {
    if (s=randomBit(p))
      if (s==1) vals[i] += r;   // add r (no need for modular reduction)
      else      SubMod(vals[i], vals[i], r, M);
    if (m>2 && (s=randomBit(p)))
      if (s==1) AddMod(vals[i], vals[i], rSqr, M); // add r^2
      else      SubMod(vals[i], vals[i], rSqr, M); // subtract r^2
  }
  if (m>4) {
    r_to_m = rSqr;
    for (j=4; j<m; j*=2) {
      SqrMod(r_to_m, r_to_m, M); // r^j := (previous-r^j)^2
      for (i=0; i<n; i++)        // vals[i] += b_{i,j} * r^j mod M
	if (s=randomBit(p))
	  if (s==1) AddMod(vals[i],vals[i],r_to_m,M);
	  else      SubMod(vals[i],vals[i],r_to_m,M);
    }
  }
  else if (m<4) { // if m==2 or 3 we're done, just return the correct r_to_m
    if (m==2) r_to_m = rSqr;
    else      MulMod(r_to_m,rSqr,r,M);
    return;
  }

  // Handle all the other powers of r

  // compute r^j,r^{2j},r^{4j},..., and add to all values
  ZZ r_odd_pwr = r;
  for (int j=3; j<m; j+=2) {
    MulMod(r_odd_pwr, r_odd_pwr, rSqr, M); // next odd power of r
    r_to_m = r_odd_pwr;
    k = j;
    while (true) {
      for (i=0; i<n; i++)     // vals[i] += b_{i,j} * r^k mod M
	if (s=randomBit(p))
	  if (s==1) AddMod(vals[i],vals[i],r_to_m,M);
	  else      SubMod(vals[i],vals[i],r_to_m,M);
      k *= 2;
      if (k >= m) break;
      SqrMod(r_to_m, r_to_m, M); // r^k := (previous-r^k)^2 mod M
    }
  }

  // r_odd_power is r^{m-1} or r^{m-2}, depending  on whether m is even or odd
  if (m&1) MulMod(r_to_m, r_odd_pwr, rSqr, M);
  else     MulMod(r_to_m, r_odd_pwr, r, M);
}
