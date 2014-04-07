// Copyright (C) IBM, All Rights Reserved

/*** randomStuff.h - utility functions
 *
 *   - random_NormalZZ(): normal distrubition integer [[not implemented]]
 *   - random_ZZ(...):    either uniform or normal, depending on flag
 *   - random_vecZZ(...): random integer vector
 *   - random_ZZX(...):   random integer polynomial
 ********************************************************************/
#ifndef _RANDOMSTUFF_H_
#define _RANDOMSTUFF_H_
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZX.h>

// Define flags
#define IS_GAUSSIAN 1
#define IS_MONIC 2

// void random_NormalZZ(ZZ &i, long logStdev); [[not implemened yet]]

inline void random_ZZ(ZZ &i, long bitsize, int flags=0) {
  // if (IS_GAUSSIAN & flags) random_NormalZZ(i,bitsize); // choose as Gaussian
  // else {
    RandomBits(i,bitsize);      // uniform in [0,2^t-1]
    i -= power2_ZZ(bitsize-1);  // uniform in [-2^{t-1}, 2^{t-1}]
  // }
}

void random_vecZZ(vec_ZZ &vec, long dim, long bitsize, int flags=0);

inline void random_ZZX(ZZX &poly, long dim, long bitsize, int flags=0)
{
  clear(poly);
  random_vecZZ(poly.rep, dim, bitsize, flags);
  if (IS_MONIC & flags) SetCoeff(poly, dim-1); // Set leading coeff to 1
  else poly.normalize();       // make sure that leading coef is nonzero
}


// returns -1/0/1, where abs(randombit()) is a Bernoulli-p bit
inline int randomBit(double p=0.5)
{
  static unsigned long f = (unsigned long) -1;
  unsigned long w = RandomWord();
  return (w <= p*f)? ((w&1)? 1: -1): 0;
}

#endif /* ifndef _RANDOMSTUFF_H_ */
