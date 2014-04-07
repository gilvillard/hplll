// Copyright (C) IBM, All Rights Reserved

#include <NTL/ZZ.h>

NTL_CLIENT
#include "fhe.h"

/** The decryption algorithm, with pubkey (d,r) seckey w and ctxt c, is:
 *
 * 1. Compute e = c*w mod d
 * 2. Return e mod 2
 ********************************************************************/
unsigned int FHEkeys::decrypt(const FHEctxt& c) const
{
  ZZ e;
  MulMod(e, c, w, det); // Compute e = c*w mod d

  // We need to map it to (-d/2,+d/2], but only care about the LSB.
  // To avoid allocating memory for another ZZ, we multiply e itself
  // by two, then check if it is bigger than det.

  unsigned int lsb = bit(e,0);
  if (IsOdd(det)) { // subtracting det would have changed the LSB
    e <<= 1;
    if (e>det) lsb ^= 1; // toggle the bit
  }
  return lsb;
}
