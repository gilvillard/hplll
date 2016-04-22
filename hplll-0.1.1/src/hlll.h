/* Householder LLL 

Created Mar 18 jan 2011 18:08:24 CET  
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


#ifndef HPLLL_HLLL_H
#define HPLLL_HLLL_H


#include  "defs.h"
#include  "mat.h"
#include  "matpe.h"
#include "matmixed.h"




namespace hplll { 

// MatrixFT pour  matrix<FP_NR<FT> > 
// MatrixZT pour  matrix<Z_NR<ZT> >

template<class ZT, class FT, class MatrixZT, class MatrixFT>
class Lattice
{

 protected:

  MatrixZT B;

  MatrixZT U;

  int n,d; 

  bool transf;

  unsigned int nblov_max;

  
  int nmaxkappa;

  // Floating point objects concerned a possible precision change 
  // ************************************************************

  MatrixFT R; 
  MatrixFT Rkept; 

  MatrixFT V; 

  MatrixFT Bfp; // floating point B 

  matrix<FP_NR<FT> > VR; // Difference between MatrixFT and matrix<FP_NR<FT> >  in Exp 
  vector<FP_NR<FT> > normB2; // Square norm  

  FP_NR<FT> x; // For size reduction 

  vector<FP_NR<FT> > toR; // Some assignment in householder_r 

  vector<int> structure; 

  vector<int> col_kept;

  vector<int> kappamin;
 
  vector<int> descendu;

public:


  // Timings 
  // ******* 
  unsigned int tps_reduce;
  unsigned int tps_householder;
  unsigned int tps_prepare;
  unsigned int tps_swap;
  unsigned int nblov,nbswaps;
  unsigned int tps_redB;
  
  int compteur;   // while counting 
  int tmpcompt;   // Debug or test counting 

  int householder_r(int kappa); 
  int householder_v(int kappa); 

  int householder();
  
  int hsizereduce(int kappa, int fromk=0);
  
  int decrease(int kappa);
  int seysenreduce(int kappa);
  int seysen_flag;

  int fast_long_flag;

  int hlll(double delta, bool verbose=false);
   
  void assignL(ZZ_mat<mpz_t> L_in);

  unsigned int setprec(unsigned int prec);
  unsigned int getprec();

  unsigned int set_nblov_max(unsigned int nb); 

  ZZ_mat<ZT> getbase();

  MatrixZT getmixedbase();

  ZZ_mat<ZT> getU();

  // Not MatrixFT for the exp case 
  matrix<FP_NR<FT> > getR(); 

  Lattice(ZZ_mat<ZT> A, bool forU=false, int reduction_method=0, int long_flag = 1); 

  Lattice(matrix<FP_NR<mpfr_t> > F, ZZ_mat<ZT> A, bool forU, int reduction_method);

  //Lattice(MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> > A, bool forU=false, int reduction_method=0);

  //Lattice(ZZ_mat<ZT> A, long t, long sigma, bool forU, int reduction_method); 

  void init(int n, int d, bool forU=false);

  void assign(ZZ_mat<ZT> A);

  void assign(MatrixZT A); 

  void shift_assign(ZZ_mat<ZT> A,  vector<int> shift, int sigma);

  void put(ZZ_mat<ZT> A, long upperdim, long t, long sigma=0); 

  void mixed_put(MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> > A, long t, long sigma=0); 

  void shift(ZZ_mat<ZT> A, long m, long sigma); 

  void shiftRT(long sigma); 

  // Only in the mpfr case (when possible to change the precision)
  // *************************************************************

  void isreduced(double delta);

#define ANY 0
#define TRIANGULAR_PROPER 1   

#define DEFAULT_PREC 0 
#define UNKNOWN_PREC 1
#define CHECK  1

// PREC IF >=2

  long lcond(int tproper =0, int flagprec=0, int flagheur=0);


  //~Lattice();
};

} // end namespace hplll


#include "hlll.cc"


#endif
