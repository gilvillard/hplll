/* Integer FGAS decomposition 

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


#ifndef HPLLL_DECOMPZ_H
#define HPLLL_DECOMPZ_H 

#include  "defs.h"
#include  "mat.h"
#include  "matpe.h"

// See links with decomp.h 

// MatrixZT pour  matrix<Z_NR<ZT> >    F, W 
// MatrixFT pour  matrix<FP_NR<FT> >   R 

template<class ZT, class FT, class MatrixZT, class MatrixFT>
class ZFgas
{
 protected:

 // Input integer matrix 
  MatrixZT F;   // Copy of the matrix d x n, is modified 

  MatrixZT U; // Transformation matrix n x n 
  MatrixZT W; // Inverse transpose transformation matrix n x n

  int d; 
  int n;   

  bool transf;

  MatrixFT R;
  MatrixFT V;

  vector<FP_NR<FT> > normB2; // Squared norm  

public:

  long int dec;  // Global variable dec used for switching between HJLS and PSLQ
                 //   dec is the dimension e.g. dec=1 one vector of the space 
                 //   for simultaneous relations detection 

  int nblov;   // Number of iterations 

  double confidence_gap;  // Gap/ratio between non zero diag and zero detection 

  int ldim;   // For the trapeze form, lattice dimension already found   

  int tmpcompt;   // Dummy counting 

  int householder();

  int householder(int kappa);  

  // Should not differ too much from the hlll size reduction function 
  int hsizereduce(int kappa);
  
  // Main decomposition function 
  int decomp(long double gamma, long int targetdim);


  // A VOIR POUR F ????????????????
  ZZ_mat<ZT> getbase();

  ZZ_mat<ZT> getU();
  ZZ_mat<ZT> getV();

  // Not MatrixFT for the exp case 
  matrix<FP_NR<FT> > getR(); 

  // Construction 
  // ------------

  // The precision can be different from the outside global one 
  // Attention !!! Put back the previous one if needed since 
  // this is done through  mpfr_set_default_prec(prec); 

  ZFgas(ZZ_mat<ZT> Finput, int forUV, long int dec);  

  void init(int d, int n, int forUV, long int dec);


  unsigned int setprec(unsigned int prec); 

  unsigned int getprec(); // "External prec, input fgas

  int set_confidence_gap(double ratio);

  unsigned int setprec_internal(long prec);
};

#endif

