/* OpenMP LLL 

Created Mer  9 avr 2014 14:07:39 CEST 
Copyright (C) 2014  Gilles Villard 

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


#ifndef HPLLL_SLLL_H
#define HPLLL_SLLL_H

#include "hlll.h"

namespace hplll { 


template<class ZT, class FT, class MatrixZT, class MatrixFT>
class PLattice
{
 protected:

  MatrixZT B;

  MatrixZT RZ;

  ZZ_mat<ZT> U;

  MatrixZT Uglob;

  bool transf;

  int n,d; 


  // Floating point objects concerned a possible precision change 
  // ************************************************************

  MatrixFT R; 

  MatrixFT Rt; 

  MatrixFT V;

 public:

  int nblov;

  int hlll(double delta, int K, int level=0, unsigned int lovmax=4294967295);
 
  void even_hsizereduce(int S, int prec, bool refresh); // Householder is refreshed or not 
  void odd_hsizereduce(int S, int prec);

  unsigned int setprec(unsigned int prec);
  unsigned int getprec();

  long  approx_cond();
  int householder();
  int phouseholder(int S, int prec);

  ZZ_mat<ZT> getU();

  ZZ_mat<ZT> getbase();

  // Not MatrixFT for the exp case 
  matrix<FP_NR<FT> > getR(); 

  SLattice(ZZ_mat<ZT> A, bool forU=false, int reduction_method=0); 

  void init(int n, int d, bool forU);

  //~SLattice();
};


} // end namespace hplll

#include "slll.cc"

#endif 

