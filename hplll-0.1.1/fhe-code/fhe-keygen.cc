// Copyright (C) IBM, All Rights Reserved

#include <cstring>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

NTL_CLIENT
#include "randomStuff.h"
#include "fhe.h"

static int InvModFx(ZZ& res, ZZ& root, ZZ& wi,
		    const ZZX& q, unsigned long logn);
static void Gz_Mod_z2(ZZ& g0, ZZ& g1, const ZZX& q, unsigned long logn);

/** The keyGen algorithm (with bitsize=t and dimension=N) is:
 *
 * 1. Choose a random vector as v=<r0 r1 ... r{N-1}> + <0 2^t 0 ... 0>
 *    where the ri's are rounded Normal with variance 2^{2t}/N log N.
 *    Let v(x) be the polynomial with these coefficients.
 *
 * 2. Run an optimized algorithm to compute the resultant of v(x)
 *    and F(x)=x^N+1, as well as one an N-th root of -1 modulo the
 *    resultant, and one odd coefficient of the scaled inverse w(x),
 *    such that w(x)*v(x) = resultant mod F(x).
 *
 *    Denote the coefficients of w(x) by <w0 w1 ... w{N-1}>, then it
 *    should hold that w1/w0 = w2/w1 =...= w{n-1}/w{n-2} = -w0/w{n-1}.
 *    The underlying public key is (resultant=D,root=R) and the secret
 *    key is one odd coefficient of w (denoted W).
 *
 * 3. Add "squashing info" to the secret key: a big set of numbers
 *    such that a sparse subset of it sums up to W modulo D, and an
 *    encryption of the charateristic vector of that subset.
 ********************************************************************/
int FHEkeys::keyGen()
{
  /* The class FHEkeys has data fields: params, det,root,w, and also *
   * the "squashing info" which is an std::vector pkBlocks.          */
  int i;
  double timer;

  // unsigned long ell = prms.logn + NextPowerOfTwo(prms.logn);// log(n log n)
  unsigned long bitsize = prms.t; // - (ell/2);
  // ZZ T = power2_ZZ(prms.t);

  unsigned long n = 1UL << prms.logn;

  int nTrials=0;
  ZZX v(INIT_SIZE, n);
  timer = -GetTime();
  do {                        // Repeat until you find good key
    nTrials++;

    // 1. Choose a random polynomial, the coefficient of the generating
    // polynomial are chosen in the interval [+- 2^t]
    random_ZZX(v, n, bitsize, IS_GAUSSIAN); // IS_GAUSSIAN not yet implemented

    //    v.rep[1] += T;

    // make sure that the coefficient sum of v is odd
    long a = trunc_long(coeff(v,0), 1); // the LSB of the coefficient
    for (int i=1; i<=deg(v); i++) {
      a ^= trunc_long(coeff(v,i), 1);
    }
    if ((a&1)==0) v.rep[1]++;

    // 2. Compute the resultant and one coefficient of the scaled inverse
  } while (!InvModFx(det,root,w,v,prms.logn));
          // InvModFx returns true on success, false on failure

#ifndef NO_TIMER
  timer += GetTime();
  cerr << "Computed det, root, w in "<< nTrials <<" trials, in time "
       << timer << "\n"; 
#endif

  /* Choose prms.s arithmetic progressions modulo det, each with
   * prms.S elements and ratio prms.R, such that there is a set
   * of one elements of each progression that sums up to w
   */
  {ZZ sum, factor;
  pkBlocks.resize(prms.s); // allocate memory for the s blocks

  // choose at random the first s-1 progressions
  sum = 0;
  for (i=0; i<prms.s-1; i++) {
    RandomBnd(pkBlocks[i].x, det); // choose a random element in [0..det-1]
    pkBlocks[i].idx = RandomBnd(prms.S);      // and an index in [0..S-1]

    // add R^{idx} * x to the sum
    power2(factor, pkBlocks[i].idx * prms.logR); // 2^{idx*logR} = R^{idx}
    if (factor>det) factor %= det;
    MulMod(factor, pkBlocks[i].x, factor, det);    // factor *= x   mod det
    AddMod(sum, sum, factor, det);                 // sum += factor mod det
  }

  // choose the last progression so as to get the sum to equal w
  i = prms.s -1;
  pkBlocks[i].idx = RandomBnd(prms.S);        // a random index in [0..S-1]
  SubMod(sum, w, sum, det);                   // sum = w - sum mod det
  power2(factor, pkBlocks[i].idx * prms.logR);// 2^{idx*logR} = R^{idx}
  if (factor>det) factor %= det;
  InvMod(factor, factor, det);                // factor = R^{-idx} mod det
  MulMod(pkBlocks[i].x, sum, factor, det);    // sum * R^{-idx} mod det

#ifndef NO_DEBUG
  // Sanity-check: verify that \sum_i x_i R^{idx_i} = w mod det
  sum = 0;
  for (i=0; i<prms.s; i++) {
    power2(factor, pkBlocks[i].idx * prms.logR); // 2^{idx*logR} = R^{idx}
    if (factor>det) factor %= det;
    MulMod(factor, pkBlocks[i].x, factor, det);    // factor *= x   mod det
    AddMod(sum, sum, factor, det);                 // sum += factor mod det    
  }
  if (sum != w) Error("FHEkeys::recrypt: sum_i x_i R^{idx_i} != w mode det");
#endif
  }
  // Now encrypt the characteristic vector of the idx'es (in a compact form)

  // Compute the number of ciphertext for each progression,
  // i.e., an integer N such that N(N-1)/2 > S (see fhe.h)
  unsigned long nCtxts = mChoose2(prms.S);

  // initialize bits to all zero, then set some of them to one
  unsigned int bits[nCtxts*prms.s];
  memset(bits, 0, sizeof(bits));
  for (i=0; i<prms.s; i++) {
    unsigned long j1, j2;
    // let j1,j2 be the idx'th pair in {nCtxts choose 2}
    encodeIndex(j1, j2, pkBlocks[i].idx, nCtxts);
    bits[i*nCtxts +j1] = bits[i*nCtxts +j2] = 1; // set these two bits to one
  }
#ifndef NO_TIMER
  cerr << (nCtxts*prms.s) << " secret-key bits ";
  timer = -GetTime();
#endif
  encrypt(ctxts, bits, nCtxts*prms.s);        // encrypt all secret key bits
#ifndef NO_TIMER
  timer += GetTime();
  cerr << "encrypted in time " << timer << "\n";
#endif
  return nTrials;
}

/********************************************************************/
/******                   Utility functions                    ******/
/********************************************************************/

/* Polynomial inversion modulo X^N+1, where N is a power of two
 *
 * This procedure returns the resultant r of q(x) and f(x)=x^N+1,
 * and also one coefficient from the scaled inverse qInv(x), such
 * that q(x)*qInv(x) = r (mod f(x)).
 */
int InvModFx(ZZ& res, ZZ& root, ZZ& wi, const ZZX& q, unsigned long n) 
{
  int i;
  unsigned long N = 1UL << n;
  // GV unused ? 
  //double invTime = 0.0, chkTime = 0.0;

  // Sanity-check
  if (deg(q)>=N) Error("MyInvMod: deg q >= deg f");

  // compute resultant and w0
  ZZ w0,w1;
  Gz_Mod_z2(res, w0, q, n);

  if (!IsOdd(res)) {   // Resultant must be odd
#ifndef NO_DEBUG
    cerr << "InvModFx: resultant is even!\n";
#endif
    return 0;
  }

  // repreat for the polynomial x*q(x) mod x^N+1
  {ZZX qx(INIT_SIZE,N);
  for (i=N-1; i>0; i--) SetCoeff(qx, i, coeff(q,i-1)); // copy 1st N-1 coeffs
  NTL::negate(qx.rep[0], coeff(q,N-1));                // negate last coeff
  qx.normalize();
  Gz_Mod_z2(res, w1, qx, n);}

  // now that we have res, w0, w1, set root = w1/w0 mod res

  // First some NTL conversions, ensuring that things are positive
  if (sign(res)==-1) {
    NTL::negate(res,res);
    NTL::negate(w0,w0);
    NTL::negate(w1,w1);
  }
  if (sign(w0)<0) w0 += res;
  if (sign(w1)<0) w1 += res;

  {ZZ tmp;
  if (InvModStatus(tmp, w0, res) != 0) {
#ifndef NO_DEBUG
    cerr << "InvModFx: w0 not invertible!\n";
#endif
    return 0; // verify that w0^{-1} exists
  }
  MulMod(root, tmp, w1, res);        // root= w1 * w0^{-1} mod res
  PowerMod(tmp, root, N, res);
  if (tmp+1 != res) {
#ifndef NO_DEBUG
    cerr << "InvModFx: root^N != -1 mod res!\n";
#endif
    return 0;
  }}

  // Let wi be any odd coefficient of W (when mapped to [-res/2,+res/2])

  if ( ((w0<=res/2)&&IsOdd(w0)) || ((w0>res/2)&&!IsOdd(w0)) )
    wi = w0;
  else if ( ((w1<=res/2)&&IsOdd(w1)) || ((w1>res/2)&&!IsOdd(w1)) )
    wi = w1;
  else for (i=2; i<N; i++) {
    MulMod(w1, w1, root, res);
    //    cerr << "w"<<i<<" = " << w1 << "\n";
    if ( ((w1<=res/2)&&IsOdd(w1)) || ((w1>res/2)&&!IsOdd(w1)) ) {
      wi = w1;
      break;
    }
  }

  return ((i==N) ? 0: 1); // We get i==N only if all the wi's are even
}

/* Compute the first two coefficients of the polynomial 
 *         G(z) = \prod_i (q(ri) - z).
 * The ri's are all the N-th roots of -1 (over the complex numbers),
 * where N = 2^n.  Returns in g0 the constant term, and in g1 the
 * linear term divided by N.
 */
static void Gz_Mod_z2(ZZ& g0, ZZ& g1, const ZZX& q, unsigned long n)
{
  int i;
  unsigned long N = 1UL << n;

  ZZX V(q);                                            // V = q
  ZZX U(INIT_SIZE, N);   set(U);                       // U = 1
  ZZX F(INIT_SIZE, N+1); SetCoeff(F,0); SetCoeff(F,N); // F(x) = x^N +1
  ZZX V2(INIT_SIZE, N);

  while (N>1) {
    V2 = V;
    for (i=1; i<=deg(V2); i+=2) {       // set V2(x) := V(-x)
      NTL::negate(V2.rep[i],V2.rep[i]); // negate odd coefficients
    }
    V2.normalize();

    MulMod(V, V, V2, F); // V := V(x) * V(-x) mod f(x)
    MulMod(U, U, V2, F); // U := U(x) * V(-x) mod f(x)

    // Sanity-check: verify that the odd coefficients in V are zero
    for (i=1; i<=deg(V); i+=2) if (!IsZero(V.rep[i])) {
      Error("MyInvMod: odd coefficient in V(x)*V(-x) not zero!");
    }
    // "Compress" the non-zero coefficients of V
    for (i=1; i<=deg(V)/2; i++) V.rep[i] = V.rep[2*i];
    for (   ; i<=deg(V);   i++) clear(V.rep[i]);
    V.normalize();

    // Set U to the "compressed" ( U(x) + U(-x) ) /2
    for (i=0; i<=deg(U)/2; i++) U.rep[i] = U.rep[2*i];
    for (   ; i<=deg(U);   i++) clear(U.rep[i]);
    U.normalize();

    // Set N := N/2 and update F accordingly
    SetCoeff(F,N,0);
    N >>= 1;
    SetCoeff(F,N);
    F.normalize();
  }

  // The free term of V is g0, that of U is g1
  g0 = ConstTerm(V);
  g1 = ConstTerm(U);
}
