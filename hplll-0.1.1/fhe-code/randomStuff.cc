// Copyright (C) IBM, All Rights Reserved

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>

NTL_CLIENT
#include "randomStuff.h"

/** A procedure that fills an integer vector with either uniform
 * signed integers with bitsize bits, or rounded Gaussians of
 * variance 2^{2*bitsize}
 ********************************************************************/
void random_vecZZ(vec_ZZ &v, long dim, long bitsize, int flags)
{
  /* random_ZZ() is defined in fhe.h, it either calls random_NormalZZ
   * or the NTL functrion RandomBits, depending on the isGaussian flag
   */
  v.SetLength(dim);
  for (int i=0; i<v.length(); i++) random_ZZ(v[i], bitsize, flags);
}
