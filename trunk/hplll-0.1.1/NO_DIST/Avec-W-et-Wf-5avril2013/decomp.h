/* FGAS decomposition 

Created Created Jeu  4 oct 2012 11:22:29 CEST 
Copyright (C) 2011, 2012, 2013      Gilles Villard 

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


#ifndef HPLLL_DECOMP_H
#define HPLLL_DECOMP_H


#include  "defs.h"
#include  "mat.h"
#include  "mat-exp.h"

// MatrixRT pour  matrix<FP_NR<FT> >   F 
// MatrixZT pour  matrix<Z_NR<ZT> >   U 
// MatrixFT pour  matrix<FP_NR<FT> >   R 

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
class Fgas
{
 protected:

 // Input floating point matrix possible concerned with the precision (and possible precision change) 
  MatrixRT F;   // Copy of the matrix d x n, is modified 

  MatrixRT F0;   // Copy of the matrix d x n, is kept 

  MatrixRT Y;   // Copy of the matrix d x 1 A REGLER, is kept  

  MatrixRT Wf; // Real number version of the transform W

  MatrixZT U; // Transpose inverse of the transformation matrix n x n 
  MatrixZT W; // Transformation matrix n x n 

  int d; 
  int n;   

  MatrixFT R;
  MatrixFT V;

  vector<FP_NR<FT> > normB2; // Squared norm  

public:

  int nblov;

  int ldim;   // For the trapeze form, lattice dimension already found   

  int tmpcompt;   // Dummy counting 

  int householder();

  int householder(int kappa);  

  // Should not differ too much from the hlll size reduction function 
  int hsizereduce(int kappa);

  
  // Main decomposition function 
  int decomp(long double gamma, long int targetdim, long int dec);


  // A VOIR POUR F ????????????????
  ZZ_mat<ZT> getbase();

  ZZ_mat<ZT> getU();
  
  // Not MatrixFT for the exp case 
  matrix<FP_NR<FT> > getR(); 

  // Construction 
  // ------------

  // The precision can be different from the outside global one 
  // Attention !!! Put back the previous one if needed since 
  // this is done through  mpfr_set_default_prec(setprec); 

  Fgas(MatrixRT Finput, long setprec);  

  void init(int d, int n, long setprec);


  // Does not do anything here (see L1)
  // ----------------------------------

  unsigned int setprec(unsigned int prec); 

};

#endif

