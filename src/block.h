/* Very preliminary implementation of recursive L1, knapsack case

Created Created   Mar 13 mar 2012 10:37:59 CEST
Copyright (C) 2012, 2013      Gilles Villard

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

#include "hlll.h"

#ifndef HPLLL_BLOCK_H
#define HPLLL_BLOCK_H

namespace hplll {

// One level : block size bs
// Knapsack cf put = 1
// B resize outside !!!
void blevel(ZZ_mat<mpz_t>& B, const ZZ_mat<mpz_t> Ain, int bs, int dec, int lllmethod = HLLL) {

  int n, d, k, i, j;


  n = Ain.getRows();
  d = Ain.getCols();

  // A block matrix
  ZZ_mat<mpz_t> A;
  A.resize(n, bs);

  // !!! For forcing dense structure !!!!
  for (i = 0; i < n; i++)
    for (j = 0; j < bs; j++)
      A(i, j) = 1;

  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L(A, NO_TRANSFORM, DEF_REDUCTION);



  // Loop over the blocks
  for (k = 0; k < d / bs; k++) {

    for (i = 0; i < n; i++)
      for (j = 0; j < bs; j++)
        A(i, j) = Ain(i, k * bs + j);


    if (k >= dec) { // We skip the first blocks
      L.assign(A);
      L.hlll(0.99);

      A = L.getbase();
    }

    for (i = 0; i < n; i++) {
      for (j = 0; j < bs; j++) {
        B(i, k * bs + j) = A(i, j);
      }
    }


  } // end over the blocks
  cout << "nblov: " << L.nblov << endl;
};




} // end namespace hplll


#endif
