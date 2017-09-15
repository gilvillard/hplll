/* LLL based relation algorithms

Created Jeu  7 mar 2013 15:01:41 CET
        Ven 27 mar 2015 14:18:38 CET

Copyright (C) 2013      Gilles Villard
Modified Mar 12 sep 2017 15:34:12 CEST

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



#ifndef HPLLL_RELATIONS_H
#define HPLLL_RELATIONS_H


#include  "hlll.h"



namespace hplll {


template<class ZT, class FT>
class FPTuple
{

protected:


  int m, d;

  vector<FP_NR<mpfr_t> > fpv;


  // Dirty templating, for matrix type in calling hplll MatrixPE or matrix, should disappear

  bool RELATIONS_DPE_FLAG = false;
  bool RELATIONS_LDPE_FLAG = false;

  // Because of restrictions on input types for fplll, should be temporary

  int call_fplll_f(ZZ_mat<FT> &b, ZZ_mat<FT> &u, double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,          \
                   LLLMethod method = LM_WRAPPER, FloatType floatType = FT_DEFAULT,               \
                   int precision = 0, int flags = LLL_DEFAULT);

  int call_fplll_z(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,          \
                   LLLMethod method = LM_WRAPPER, FloatType floatType = FT_DEFAULT,               \
                   int precision = 0, int flags = LLL_DEFAULT);

public:


  FPTuple(vector<FP_NR<mpfr_t> > fpvin);

  //~FPTuple();


  int relation_f(ZZ_mat<mpz_t>& C, long alpha,
                 long confidence_gap = 60, long shift = 200, long increment = 20,
                 int lllmethod = FPLLL, double delta = 0.99);


  int relation_f_z(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A,  int alpha,
                   long confidence_gap = 60, long shift = 200, long increment = 20,
                   int lllmethod = FPLLL, double delta = 0.99);

  int detect_lift_f_z(ZZ_mat<ZT>& U, ZZ_mat<mpz_t> L_in, ZZ_mat<FT> A_in_f, int& new_def, int def,
                      int target_def, FP_NR<mpfr_t>& new_quot,
                      long confidence_gap = 60, long shift = 200, long increment = 20,
                      int lllmethod = FPLLL, double delta = 0.99);



  int relation_z(ZZ_mat<mpz_t>& C, long alpha,
                 long confidence_gap = 60, long shift = 10, int truncate = -1, int lllmethod = FPLLL, double delta = 0.99);



  int relation_lll_z(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A, long alpha,
                     long confidence_gap = 60, long shift = 10, int truncate = -1, int lllmethod = FPLLL, double delta = 0.99);




};

} // end namespace hplll


#include "relations.cc"


#endif
