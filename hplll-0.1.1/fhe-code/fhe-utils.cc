// Copyright (C) IBM, All Rights Reserved

#include <iostream>
#include <NTL/ZZ.h>

NTL_CLIENT
#include "fhe.h"

/********************************************************************
 ******************       Utility functions        ******************
 ********************************************************************/

// returns in j1,j2 the i'th pair (in lexicographic order) in {m choose 2}.
// Indexing starts at zero, so we have m-1 pairs of the form (0,j2),
// m-2 pairs of the form (1,j2), etc. The index i is also zero-based,
// so i=0,1,...m-2 is encoded as (0,1),(0,2),...,(0,m-1), then
// i=m-1,m,m+1,...,2m-4 is encoded as (1,2),(1,3),(1,4),...,(1,m-1), etc.
void FHEkeys::encodeIndex(unsigned long& j1, unsigned long& j2,
			  unsigned long i, unsigned long m)
{
  if (i > m*(m-1)/2) { // sanity check
    cerr << "encodeIndex: index "<< i
	 << " can't be encoded in ("<<m<<" choose 2)\n";
    abort();
  }
  for (j1=0; j1<m-1; j1++) {
    unsigned long pairsFor_j1 = m-j1-1;
    if (pairsFor_j1 > i) { // i is one of the pairs (j1,*)
      j2 = i + j1+1;       // j2 is in j1+1,...,m-1
      return;
    }
    else i -= pairsFor_j1; // we already counted the pairs for this j1
  }
  // we should never get to this point
  NTL::Error("encodeIndex: something went wrong, maybe an overflow?");
}

// Verify that c*w is close enough to a multiple of det
long FHEkeys::testCtxt(const FHEctxt& c) const
{
  ZZ e;
  MulMod(e, c, w, det);        // Compute e = c*w mod d
  if (e > det/2) sub(e,det,e); // e = abs(distance from e to multiple of det)
  return NumBits(det)-NumBits(e);
}

/********************************************************************
 ******************         I/O functions          ******************
 ********************************************************************/
std::istream& operator>>(std::istream& s, FHEkeys& keys)
{
  keys.inputKeys(s);
  return s;
}

void FHEkeys::inputKeys(std::istream& s, int inIdxs, int inW)
{
  s >> prms.secprm;  // get the parameters
  s >> prms.mu;
  s >> prms.s;
  s >> prms.S;
  s >> prms.p;
  s >> prms.t;
  s >> prms.logn;
  s >> prms.logR;
  s >> prms.noise;

  // The det, root, and w of the underlying scheme
  s >> det;
  s >> root;
  if (inW) s >> w;

  // The arithmetic progressions
  pkBlocks.resize(prms.s);
  for (int i=0; i<pkBlocks.size(); i++) { 
    s >> pkBlocks[i].x;
    if (inIdxs) s >> pkBlocks[i].idx;
  }
  // finally the encrypted secret-key bits
  if (inIdxs || !inW) s >> ctxts;
}

std::ostream& operator<<(std::ostream& s, FHEkeys& keys)
{
  keys.outputKeys(s);
  return s;
}

void FHEkeys::outputKeys(std::ostream& s, int outIdxs, int outW)
{
  s << prms.secprm;  // output all the parameters
  s << " " << prms.mu;
  s << " " << prms.s;
  s << " " << prms.S;
  s << " " << prms.p;
  s << " " << prms.t;
  s << " " << prms.logn;
  s << " " << prms.logR;
  s << " " << prms.noise;

  // The det, root, and w of the underlying scheme
  s << "\n" << det;
  s << "\n" << root;
  if (outW) s << "\n" << w;

  // The arithmetic progressions
  for (int i=0; i<pkBlocks.size(); i++) { 
    s << "\n" << pkBlocks[i].x;
    if (outIdxs) s << " " << pkBlocks[i].idx;
  }
  // finally the encrypted secret-key bits
  if (outIdxs || !outW) s << "\n" << ctxts;
}

void FHEparams::print(std::ostream& st) const
{
  st << "secparam = "<<secprm<<", mu = "<<mu<<"\n";
  st << "sparse-subset size s = "<<s<<", precision p = "<<p<<"\n";
  st << "big-set size S = "<<S<<", ratio R = 2^"<<logR<<"\n";
  st << "integer-size t = "<<t<<", dimension n = 2^"<<logn<<"\n";
  st << "noise parameter = " << noise << endl;
}

/********************************************************************
 ***************  Functions to handle the parameters  ***************
 ********************************************************************/
void FHEparams::setPrms(unsigned long sec, unsigned long tt, unsigned long ln)
{
  static double log2 = log(2.0);
  unsigned long deg; // degree of the squashed decryption polynomial
  unsigned long n;

  if (sec <=0) sec = 64;
  secprm=sec;
  double logsec = log((double)sec)/log2;

  // if dimenssion not specified, compute it later from default mu value
  if (ln <= 0) mu = 0.11; 

  // Compute the parameters s (sparse subset size), p (precision), d (degree)
  {unsigned long ls  // find the first power of two larger than sec/log(sec)
     = NextPowerOfTwo((unsigned long)ceil(sec/logsec));
  s = (1UL<<ls) -1;  // s is one less than that power of two
  deg = 2*s;         // deg = 2s
  p = ls;}           // p = log(s+1)

  // The big-set size S (first constraint)
  double logS = sec/(double)((s+1)/2); // S^{ceil{s/2}} >= 2^{secparam}
  S = (unsigned long) ceil(exp(logS*log2));

  // If not specified, compute t using the formulas
  if (tt > 0) t = tt;
  else {
    double t2 = ceil(deg*((logS/2) +logsec +4)) +p+1;
    // add log(sqrt(n)) to the bit-size
    if (ln > 0) t2 += (ln+1)/2; 
    else t2 += log(ceil(t2*sec/(mu*logsec)))/(2*log2) +0.5;// estimate of same
    t = (unsigned long) ceil(t2);
  }

  if (ln>0) { // If n is specified, compute mu from n,t
    logn = ln;
    n = 1UL<<logn;
    mu = t*sec/(n*logsec);
  } else {    // Else compute (log of) n from t, mu
    n = (unsigned long) ceil(t*sec/(mu*logsec));
    logn = NextPowerOfTwo(n);
    n = 1UL<<logn;
  }

  double S2 = ceil( sqrt(t*n*sec/(mu*logsec)) / s );
  if (S < S2) {
    S = (unsigned long) S2;
    // Update t,n if needed
    if (tt <= 0) {
      logS = log(S2)/log2;
      t = (unsigned long) ceil(deg*((logS/2) +logsec +4) +(logn/2) +p+1);
      if (ln <= 0) {
	n = (unsigned long) ceil(t*sec/(mu*logsec));
	logn = NextPowerOfTwo(n);
	n = 1UL<<logn;
      }
    }
  }

  // Compute (log of) the ratio R between big-set elements
  logR = (unsigned long) ceil(n*t / (double)(s*(S-1)));

  // Compute the noise level, such that {n choose noise/2} > 2^secparam
  noise = (unsigned long) ceil(2*secprm/(double)(logn-2));
  if (noise < s) noise = s;      // make sure we're not over-optimistic
}
