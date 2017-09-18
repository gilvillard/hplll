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



#ifndef HPLLL_RELATIONS_Z_CC
#define HPLLL_RELATIONS_Z_CC



namespace hplll {


/************************************************************************************

   Temporary wrapping the call to fplll for no compilation error

   e.g. until  Z_NR<__int128_t> available

**************************************************************************************/



template<> int
FPTuple<long, double, matrix<FP_NR<double> > >::call_fplll(ZZ_mat<long> &b, ZZ_mat<long> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(b, u, delta, eta, method, floatType, precision);

  return status;
}


template<> int
FPTuple<long, long double, matrix<FP_NR<long double> > >::call_fplll(ZZ_mat<long> &b, ZZ_mat<long> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(b, u, delta, eta, method, floatType, precision);

  return status;


}



template<> int
FPTuple<__int128_t, double, matrix<FP_NR<double> > >::call_fplll(ZZ_mat<__int128_t> &b, ZZ_mat<__int128_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  cerr << endl << "** Error in relations: no fplll with Z_NR<__int128_t> **" << endl;

  exit(EXIT_FAILURE);

  return 0;
}

template<> int
FPTuple<__int128_t, long double, matrix<FP_NR<long double> > >::call_fplll(ZZ_mat<__int128_t> &b, ZZ_mat<__int128_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  cerr << endl << "** Error in relations: no fplll with Z_NR<__int128_t> **" << endl;

  exit(EXIT_FAILURE);

  return 0;
}



template<> int
FPTuple<mpz_t, double, matrix<FP_NR<double> > >::call_fplll(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(b, u, delta, eta, method, floatType, precision);

  return status;

}

template<> int
FPTuple<mpz_t, long double, matrix<FP_NR<long double> > >::call_fplll(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(b, u, delta, eta, method, floatType, precision);

  return status;

}



template<> int
FPTuple<mpz_t, dpe_t, MatrixPE<double, dpe_t> >::call_fplll(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(b, u, delta, eta, method, floatType, precision);

  return status;

}


template<> int
FPTuple<mpz_t, ldpe_t, MatrixPE<long double, ldpe_t> >::call_fplll(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;

  status = lll_reduction(b, u, delta, eta, method, floatType, precision);

  return status;

}




/***********************************************************************************

   Construction

   Fix the precision inside or let outside ?

**************************************************************************************/

// template<class ZT, class FT>
// FPTuple<ZT, FT>::FPTuple(vector<FP_NR<mpfr_t> > fpvin) {

//  d = fpvin.size();

//  fpv.resize(d);

//  for (int i = 0; i < d; i++)
//    fpv[i] = fpvin[i];


// }



template<class ZT, class FT, class MatrixFT>
FPTuple<ZT, FT, MatrixFT>::FPTuple(vector<FP_NR<mpfr_t> > fpvin) {

  d = fpvin.size();

  fpv.resize(d);

  for (int i = 0; i < d; i++)
    fpv[i] = fpvin[i];


}


/***********************************************************************************

    Relation_lll

    Calls LLL with successive lifts on mpz_t and FT bases

    alpha correct bits

 TODO: tune the parameters according to input data types

  truncate = -1 : no truncation
           = 0 : automatic choice
           > 0 : this number of bits

**************************************************************************************/


template<class ZT, class FT, class MatrixFT>  int
FPTuple<ZT, FT, MatrixFT>::relation(ZZ_mat<mpz_t>& C,  long alpha,
                                    long confidence_gap,  long shift, int truncate, int lllmethod, double delta) {


  ZZ_mat<mpz_t> L;
  L.resize(1, d);

  FP_NR<mpfr_t> t;

  for (int j = 0; j < d; j++) {
    t.mul_2si( fpv[j], alpha);
    L(0, j).set_f(t);
  }

  int found;

  found = relation_lll(C, L, alpha, confidence_gap, shift, truncate, lllmethod, delta);

  return found;

}




/***********************************************************************************

    Companion to relation_lll
    Relation from an integer matrix

    Calls LLL with successive lifts on mpz_t and FT bases

    alpha correct bits

**************************************************************************************/

template<class ZT, class FT, class MatrixFT>  int
FPTuple<ZT, FT, MatrixFT>::relation_lll(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A, long alpha,
                                        long confidence_gap, long shift, int truncate, int lllmethod, double delta) {

  Timer time;

  Timer tlll;
  Timer tprod;
  tlll.clear();
  tprod.clear();

  int m, d;
  int i, j;

  m = A.get_rows();
  d = A.get_cols();

  ZZ_mat<mpz_t> A_in;
  A_in.resize(m + d, d);


  // **** m=1 for the moment

  for (j = 0; j < d; j++)
    A_in(0, j) = A(0, j);

  for (i = 0; i < d; i++)
    A_in(m + i, i) = 1;


  int bitsize = maxbitsize(A, 0, m, d);

  // For assigning the truncated basis at each step



  int def = -bitsize;

  int target_def = -bitsize + alpha;

  int found = 0;

  FP_NR<mpfr_t> quot, new_quot, tf;
  new_quot = 1.0;  // new_quot should be bigger after the first iteration of the loop
  new_quot.mul_2si(new_quot, alpha);

  FP_NR<mpfr_t> gap;
  gap = 1.0;

  FP_NR<mpfr_t> confidence;
  // For testing 1/gap < confidence
  confidence = 1.0;
  // relié, plus petit,  au shift sur S (ex 80)
  confidence.mul_2si(confidence, -confidence_gap - shift); // La baisse absolue est plus petite que le shift

  FP_NR<mpfr_t> epsilon;
  epsilon = 10.0;


  Z_NR<mpz_t> tz, maxcol;

  ZZ_mat<ZT> T, TT;

  T.resize(m + d, d);
  TT.resize(d, m + d);


  ZZ_mat<ZT> U, UT;

  U.resize(d, d);
  UT.resize(d, d);


  Lattice<ZT, FT, matrix<Z_NR<ZT> >,  MatrixFT > Bp(T, TRANSFORM, DEF_REDUCTION);



  // Main loop on the shifts
  // -----------------------

  while (def < target_def) {


    HPLLL_INFO("Current default: ", def);

    if ((target_def - def) <= shift)
      def = target_def;
    else def += shift;


    if (truncate == -1)
      lift_truncate(T, A_in, def, 0);
    else if (truncate == 0)
      lift_truncate(T, A_in, def, shift + 2 * d);
    else
      lift_truncate(T, A_in, def, truncate);


    cout << "Size T : " << maxbitsize(T, 1, d, d) << endl;


    // HLLL and DOUBLES

    if (lllmethod == HLLL) {


      Bp.assign(T);

      Bp.hlll(delta);

      matprod_in_int(A_in, Bp.getU());
      //avec long: matprod_in_si(A_in,U);
      cout << "sizeof U: " << maxbitsize(Bp.getU(), 0, d, d) << endl;


    } // end HLLL


    // FPLLL

    else if (lllmethod == FPLLL) {

      transpose(TT, T);

      setId(UT);

      time.start();
      call_fplll(TT, UT, delta, 0.51, LM_FAST, FT_DEFAULT, 0);
      time.stop();

      tlll += time;

      time.start();
      transpose(U, UT);

      matprod_in_int(A_in, U);
      //matprod_in_si(A_in,U);
      //avec long: matprod_in_si(A_in,U);
      time.stop();

      tprod += time;


    } // end FPLLL



    // Test
    // ----

    quot = new_quot;

    if (A_in(0, 0).sgn() == 0) { // For making the gap pertinent even if 0
      tz = 1;
    }
    else
      tz.abs(A_in(0, 0));
    new_quot.set_z(tz);


    maxcol.abs(A_in(0, 0));
    maxcol.mul_2si(maxcol, def);

    // cout << "L: " << A_in(0,0) << endl;
    // cout << "Af: " << maxcol << endl;

    for (i = 0; i < d; i++) {
      tz.abs(A_in(m + i, 0));
      if (tz.cmp(maxcol) == 1) maxcol = tz;
    }

    tf.set_z(maxcol);
    new_quot.div(new_quot, tf);


    gap.div(new_quot, quot);
    gap.abs(gap);


    // ICI
    //cout << "Gap : " << gap << endl;
    //cout << "Newquot : " << new_quot << endl;
    //cout << "Digits : " << def +alpha << endl;

    //if ((gap.cmp(confidence) == -1) && (new_quot.cmp(epsilon) == -1)) {
    if (gap.cmp(confidence) == -1) {
      C.resize(d, 1);
      for (j = 0; j < d; j++)
        C(j, 0) = A_in(m + j, 0);

      gap.mul_2si(gap, shift);

      cout << endl;
      HPLLL_INFO(def + bitsize, " bits used");
      HPLLL_INFO("Candidate relation found with bit confidence: ", -gap.exponent());

      cout << "LLL : " << tlll << endl;
      cout << "Products : " << tprod << endl;
      return 1;

    }


  } // End while on the shift

  HPLLL_INFO(alpha, " digits used");

// Relation bound
// --------------

  unsigned oldprec;
  oldprec = mpfr_get_default_prec();

  mpfr_set_default_prec(2 * d);

  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A_in, NO_TRANSFORM, DEF_REDUCTION);

  B.householder();

  matrix<FP_NR<mpfr_t> > R;

  R = B.getR();

  FP_NR<mpfr_t> minr, tr;

  minr.abs(R(0, 0));
  for (i = 1; i < d; i++) {
    tr.abs(R(i, i));
    if (minr.cmp(tr) > 0) minr = tr;
  }

  cerr << endl << "** No relation found with min Rii = " << minr << endl;

  mpfr_set_default_prec(oldprec);


  return found; // 0 here


};






} // end namespace hplll


#endif
