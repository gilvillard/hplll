// Copyright (C) IBM, All Rights Reserved

#include <vector>
#include <stack>
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>

NTL_CLIENT
#include "fhe.h"

// Utility functions to compute carry bits and add numbers in binay
typedef std::stack<ZZ> ZZstack;
static void evalSymPolys(vec_ZZ& out, ZZstack& vars, long deg, const ZZ& M);
static void gradeSchoolAdd(FHEctxt& c, const mat_ZZ& vars, const ZZ& M);


#ifndef NO_DEBUG
#include <NTL/RR.h>
static ZZ fhe_recrypt_sum;
static const FHEkeys *fhe_recrypt_curKey;
static bool verifyVector(const vec_ZZ& vars,
			 const FHEctxt& c, const PKblock& block);
#endif


// Perform homomorphic decryption of the cipehrtext c
void FHEkeys::recrypt(FHEctxt& c) const
{
#ifndef NO_DEBUG
  clear(fhe_recrypt_sum);
  fhe_recrypt_curKey = this;

  // Sanity-check: verify that c*w is close enough to a multiple of det
  if (testCtxt(c)<prms.p +1)
    cerr << "FHEkeys::recrypt: c * w larger than det / 2^{p+1}!!\n";
#endif

  // From each public-key block we compute an encrypted (p+1)-bit vector
  mat_ZZ vars(INIT_SIZE, prms.s, prms.p +1);
#ifndef NO_TIMER
  double timer = -GetTime();
#endif
  for (long i=0; i<prms.s; i++) processBlock(vars[i], c, i);
#ifndef NO_TIMER
  timer += GetTime();
  cerr << "Processed "<<prms.s<<" public-key blocks in time "<<timer<<endl;
#endif

#ifndef NO_DEBUG
  ZZ tmp;
  MulMod(tmp, c, w, det);
  if (tmp != fhe_recrypt_sum)
    Error("FHEkeys::recrypt: sum from blocks differ from w*c");
#endif

  // Use the grade-school algorithm to add up these s (p+1)-bit numbers,
  // return in c the encryption of the XOR of the two left bits
#ifndef NO_TIMER
  timer = -GetTime();
#endif
  gradeSchoolAdd(c, vars, det);
#ifndef NO_TIMER
  timer += GetTime();
  cerr << "Grade-school addition algorithm in time "<<timer<<endl;
#endif
}


/* Takes two integer 0<=n<d, and returns an nBits-bit integer, with
 * the top bit viewed as the bit to the left of the binary point,
 * and the lower bits are viewed as being to the right of the point.
 * 
 * The returned value is the binary representation of the rational
 * number n/d, rounded to precision nBits-1.
 ********************************************************************/
static unsigned long getBinaryRep(ZZ n, const ZZ& d, long nBits)
// It is assumed that nBits fit in one long integer.
{
#ifndef NO_DEBUG
  RR ratio = to_RR(n)/to_RR(d);
#endif
  // (n * 2^nBits)/d gives nBits bits of precision (one more than needed)
  n <<= nBits;
  n /= d;         // integer division implies truncation

  unsigned long sn = to_long(n); // a single precision variant
  sn = (sn >> 1) + (sn & 1);     // one less bit, but round instead of truncate
     // NOTE: the addition of (sn&1) could make sn as large as 2^{nBits}. For
     // this case we need to remember the bit to the left of the binary point.

#ifndef NO_DEBUG
  // Sanity-check: check the distance between sn/2^nBits and the
  // real number n/d, verify that it is at most 1/2^{nBits}
  long twoToP = 1L<<(nBits-1);
  double approx = ((double)sn)/((double)twoToP);

  RR error = abs(ratio-approx);  // check that |ratio-approx| is small enough
  if (error*(twoToP*2) > 1) {
    cerr << "n/d = " << ratio << endl;
    cerr << "approx n/d = " << approx << endl;
    cerr << "error is 1/"<<(1/error)<<endl;
    Error("getBinaryRep: sn/2^p too far from n/d");
  }
#endif
  return sn;
}


// Processing the i'th public-key block:
// Compute encryption of top p bits of \sum_j sigma_j*(c*x*R^j mod det)/det
// + encryption of top bit (where the LSBs are XORed in)
void FHEkeys::processBlock(vec_ZZ& vars, const FHEctxt& c, long i) const
{
#ifndef NO_TIMER
  double ptimer = 0.0;
  double atimer = 0.0;
  double mtimer = 0.0;
#endif
  long nCtxts = mChoose2(prms.S);
  unsigned long baseIdx = i * nCtxts;

  int k;
  for (k=0; k<vars.length(); k++) clear(vars[k]); // initialize to zero

  unsigned long j, j1, j2;
  ZZ factor = pkBlocks[i].x;
  MulMod(factor, factor, c, det);
  vec_ZZ psums(INIT_SIZE, vars.length()); // keep partial sums
  for (j=j1=0; j1<nCtxts-1; j1++) {       // sk-bits indexed by (j1,*) pairs
    for (k=0; k<psums.length(); k++) clear(psums[k]);  // initialize to zero

    for (j2=j1+1; j2<nCtxts; j2++) {
      // get the top bits of factor/det. The code below assumes
      // that p+1 bits can fit in one unsigned long
#ifndef NO_TIMER
      ptimer -= GetTime();
#endif
      unsigned long binary = getBinaryRep(factor, det, vars.length());
      if (IsOdd(factor))     // "xor" the LSB to column 0
	binary ^= (1UL << prms.p);
#ifndef NO_TIMER
      ptimer += GetTime();
      atimer -= GetTime();
#endif

      // For every 1 bit, add the current ciphertext to the partial sums
      for (k=0; k<psums.length(); k++) if (bit(binary, k) == 1) {
	long k2 = psums.length() -k-1;
	AddMod(psums[k2], psums[k2], ctxts[baseIdx+j2], det);
      }
#ifndef NO_TIMER
      atimer += GetTime();
#endif

      j++;              // done with this element
      if (j < prms.S) { // compute next element = current * R mod det
#ifndef NO_TIMER
	ptimer -= GetTime();
#endif
	factor <<= prms.logR;
	factor %= det;
#ifndef NO_TIMER
	ptimer += GetTime();
#endif
      }
      else break;       // don't add more than S elements
    }

    // multiply partial sums by ctxts[j1], then add to sum
#ifndef NO_TIMER
    mtimer -= GetTime();
#endif
    for (k=0; k<vars.length(); k++) {
      MulMod(psums[k], psums[k], ctxts[baseIdx+j1], det);
      AddMod(vars[k],vars[k],psums[k], det);
    }
#ifndef NO_TIMER
    mtimer += GetTime();
#endif
    if (j >= prms.S) break;
  }
  // Sanity-check: j should be at least S, else we've missed some terms
  if (j < prms.S) Error("FHEkeys::processBlock: loop finished with j<S");

#ifndef NO_TIMER
  cerr << "Computing one arithemtic progression takes time "<<ptimer<<endl;
  cerr << "Computing additions takes time "<<atimer<<endl;
  cerr << "Computing multiplications takes time "<<mtimer<<endl;
#endif

#ifndef NO_DEBUG
  // Sanity-check: verify that vars[0..,p-1] are encryptions of the
  // top p bits of the rational number (c x_i R^{idx_i} mod det)/det,
  // and vars[p] encrypts the LSB of (c x_i R^{idx_i} mod det)
  if (!verifyVector(vars, c, pkBlocks[i]))
    Error("FHEkeys::processBlock: decrypted vector does not match");
#endif
}


// Use the grade-school algorithm to add up these s (p+1)-bit numbers.
static void gradeSchoolAdd(FHEctxt& c, const mat_ZZ& vars, const ZZ& M)
{
  long i,j;
#ifndef NO_DEBUG
  // count how many 1 bits in each column of vars
  vec_long zCols(INIT_SIZE, vars.NumCols());
  for (j=0; j<vars.NumCols(); j++) for (i=0; i<vars.NumRows(); i++)
    zCols[j] += fhe_recrypt_curKey->decrypt(vars[i][j]);
  cerr << "  Column-weight in input to gradeSchoolAdd: "<<zCols<<endl;
#endif

  // Below it is more convenient to have each column of the matrix in
  // a separate stack (since we would want to push carry bits on top)
  std::vector<ZZstack> columns(vars.NumCols());
  for (j=0; j<vars.NumCols(); j++) for (i=vars.NumRows()-1; i>=0; i--)
    columns[j].push(vars[i][j]);

  vec_ZZ sp; // space to store symmetric polynomials

  // add columns from right to left, upto column -1
  for (j=vars.NumCols()-1; j>0; j--) { 
    long s = columns[j].size();
    long log = NextPowerOfTwo(s); // (log of) # of carry bits to compute
    if (log > j) log = j;     // no more carry than what can reach col 0
    if ((1L<<log) > s) log--; // no more carry than what s bits can produce

    evalSymPolys(sp, columns[j], 1L<<log, M); // evaluate symmetric polys

    // The carry bits from this column are sp[2],sp[4],sp[8]... The result
    // for that column is in sp[1] (but for most columns we don't need it)
#ifndef NO_DEBUG
    zCols[j] = fhe_recrypt_curKey->decrypt(sp[1]); // record result
#endif
    long k = 2;
    for (long j2=j-1; j2>=0 && k<sp.length(); j2--) {
      columns[j2].push(sp[k]);   // push carry bits on top of their column
      k <<= 1;
    }
  }

  // The result from column -1 is in sp[1], add to it all the bit in column 0
  c = columns[0].top();
  columns[0].pop();
  while (!columns[0].empty()) {
    AddMod(c, c, columns[0].top(), M);
    columns[0].pop();
  }
#ifndef NO_DEBUG
  zCols[0] = fhe_recrypt_curKey->decrypt(c); // record result
  cerr << "  output from gradeSchoolAdd: "<<zCols<<endl;    
#endif
  AddMod(c, c, sp[1], M);
}

static void evalSymPolys(vec_ZZ &out, ZZstack& vars, long deg, const ZZ& M)
{
  long i, j;

  out.SetLength(deg+1); // initializa to [1,0,...,0,0]
  set(out[0]);
  for (i=1; i<=deg; i++) clear(out[i]);

  ZZ tmp;
  for (i=1; !vars.empty(); i++) {  // process the next variable, i=1,2,...
    for (j=min(i,deg); j>0; j--) { // compute the j'th elem. sym. poly
      MulMod(tmp, out[j-1], vars.top(), M);
      AddMod(out[j], out[j], tmp, M); // out[j] += out[j-1] * vars.top() mod M

    // At the end of the inner loop, out[j] holds the
    // j'th symmetric polynomial in the first i variables
    }
    vars.pop();  // done with this variable
  }
}

/********************************************************************/
/********* some sanity-check functions, useful for debugging ********/
/********************************************************************/

#ifndef NO_DEBUG
// Compute the top bits of (c*x*R^{idx} mod det)/det and compare to v[]
static bool verifyVector(const vec_ZZ& vars,
			 const FHEctxt& c, const PKblock& block)
{
  unsigned int v[vars.length()];
  fhe_recrypt_curKey->decrypt(v, vars);

  ZZ factor;
  long logR = fhe_recrypt_curKey->getPrms().logR;
  long p = fhe_recrypt_curKey->getPrms().p;
  const ZZ& det = fhe_recrypt_curKey->getDet();

  power2(factor, block.idx * logR);     // 2^{idx*logR} = R^{idx}
  if (factor>det) factor %= det;
  MulMod(factor, block.x, factor, det); // factor *= x mod det
  MulMod(factor, c, factor, det);       // factor *= c mod det

  AddMod(fhe_recrypt_sum, fhe_recrypt_sum, factor, det);
  unsigned long binary = getBinaryRep(factor, det, p+1);
  if (IsOdd(factor)) v[0] ^= 1; // xor the LSB

  for (long j=0; j<=p; j++) if (bit(binary,j) != v[p-j]) {
    cerr << "binary = "<<binary<<" but v[0-4]= ["
	 <<v[0]<<","<<v[1]<<","<<v[2]<<","<<v[3]<<","<<v[4]<<"]\n";
    return false;
  }
  cerr << "verifyVector: binary = "<<binary<<" verified\n";

  return true;
}
#endif /* #ifndef NO_DEBUG */
