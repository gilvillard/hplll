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

  long emax, emin;

  emin = fpv[0].exponent();
  emax = fpv[0].exponent();


  for (int i = 0; i < d; i++) {

    emin = min(emin, fpv[i].exponent());
    emax = max(emax, fpv[i].exponent());

  }

  inputgap = emax - emin;


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
                                    long confidence_gap,  long shift, int truncate, int lllmethod,  int sizemethod, double delta) {


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
                                        long confidence_gap, long shift, int truncate, int lllmethod,  int sizemethod, double delta) {

  Timer time;

  Timer tlll;
  Timer tprod;
  Timer ttrunc;
  tlll.clear();
  tprod.clear();

  Timer tphase;

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


  Z_NR<mpz_t> tz, t, maxcol;

  long foundcol;

  ZZ_mat<ZT> T, TT;

  T.resize(m + d, d);
  TT.resize(d, m + d);


  ZZ_mat<ZT> U, UT;

  U.resize(d, d);
  UT.resize(d, d);


  Lattice<ZT, FT, matrix<Z_NR<ZT> >,  MatrixFT > Bp(T, TRANSFORM, sizemethod);


//OMP
  int S = 8;

  double en, st;
  double lllt = 0.0;
  double prodt = 0.0;


#ifdef _OPENMP
  omp_set_num_threads(S);
#endif


  // Main loop on the shifts
  // -----------------------

  while (def < target_def) {

    tphase.start();

    HPLLL_INFO("Current default: ", def);

    if ((target_def - def) <= shift)
      def = target_def;
    else def += shift;


    time.start();
    if (truncate == -1)
      lift_truncate(T, A_in, def, 0);
    else if (truncate == 0)
      lift_truncate(T, A_in, def, shift + 2 * d);
    else
      lift_truncate(T, A_in, def, truncate);
    time.stop();
    ttrunc += time;


    //  DBG  
    // if (def == -135657) {

    //   ZZ_mat<ZT> AT;

    //   AT.resize(d, m + d);

    //   transpose(AT, T);

    //   cout << AT << endl;


    // }
  


    // HLLL and DOUBLES

    if (lllmethod == HLLL) {


      time.start();
      Bp.assign(T);

      Bp.hlll(delta);
      time.stop();

      tlll += time;


      //matprod_in_int(A_in, Bp.getU());


#ifdef _OPENMP

      st = omp_get_wtime();
      pmatprod_in_int(A_in, Bp.getU(), S);
      en = omp_get_wtime();
      prodt += (en - st);

#else

      st = omp_get_wtime();
      matprod_in_int(A_in, Bp.getU());
      en = omp_get_wtime();
      prodt += (en - st);

#endif


      tprod += time;

      //avec long: matprod_in_si(A_in,U);
      cout << "sizeof U: " << maxbitsize(Bp.getU(), 0, d, d) << endl;
      cout << "sizeof basis: " << maxbitsize(A_in, 1, d + 1, d) << endl << endl;



    } // end HLLL


    // FPLLL

    else if (lllmethod == FPLLL) {

      transpose(TT, T);

      setId(UT);

      //time.start();
      st = omp_get_wtime();
      call_fplll(TT, UT, delta, 0.51, LM_FAST, FT_DOUBLE, 0);
      en = omp_get_wtime();

      lllt += (en - st);

      //time.start();
      transpose(U, UT);


#ifdef _OPENMP

      st = omp_get_wtime();
      pmatprod_in_int(A_in, U, S);
      en = omp_get_wtime();
      prodt += (en - st);

#else

      st = omp_get_wtime();
      matprod_in_int(A_in, U);
      en = omp_get_wtime();
      prodt += (en - st);

#endif

      time.stop();

      tprod += time;

      cout << "sizeof U: " << maxbitsize(Bp.getU(), 0, d, d) << endl;
      cout << "sizeof basis: " << maxbitsize(A_in, 1, d + 1, d) << endl << endl;


    } // end FPLLL



    // Test
    // ----

    /* GV Mer 20 sep 2017 12:57:27 CEST
        Mofified, the small entry may appear in any column, if not in the first one and
         if integer in fixed precision (ex long) then the lift truncate may lead to an error
         for the detection (bad truncation and crash)
    */

    quot = new_quot;

    foundcol = 0;

    if (A_in(0, 0).sgn() == 0) { // For making the gap pertinent even if 0
      tz = 1;
    }
    else
      tz.abs(A_in(0, 0));


    for (int i = 1; i < d; i++) {

      t.abs(A_in(0, i));

      if (t.cmp(tz) == -1) {
        tz = t;
        foundcol = i;
      }

    }

    new_quot.set_z(tz);


    maxcol.abs(A_in(0, foundcol));
    maxcol.mul_2si(maxcol, def);


    for (i = 0; i < d; i++) {
      tz.abs(A_in(m + i, foundcol));
      if (tz.cmp(maxcol) == 1) maxcol = tz;
    }

    tf.set_z(maxcol);
    new_quot.div(new_quot, tf);


    gap.div(new_quot, quot);
    gap.abs(gap);


    /* GV Mer 20 sep 2017 12:57:27 CEST
        added the test using the input gap otherwise the computation may terminate
        without having shifted for discovering all non zero entries
    */

    if (def > inputgap - bitsize) {

      //if ((gap.cmp(confidence) == -1) && (new_quot.cmp(epsilon) == -1)) {
      // Si epsilon mettre à la valeur max quotient des nombres au départ

      if (gap.cmp(confidence) == -1) {
        C.resize(d, 1);
        for (j = 0; j < d; j++)
          C(j, 0) = A_in(m + j, foundcol);

        gap.mul_2si(gap, shift);

        cout << endl;
        HPLLL_INFO(def + bitsize, " bits used");
        HPLLL_INFO("Candidate relation found with bit confidence: ", -gap.exponent());

        //cout << endl << "Time lll: " << tlll << endl;
        //cout << "Time products: " << tprod << endl << endl;

        cout << endl << "Time trunc: " << ttrunc << endl;
        cout << endl << "Time lll: " << lllt << endl;
        cout << "Time products: " << prodt << endl << endl;

        return 1;

      }

    }

    tphase.stop();
    cout << "--- Phase " << def << " time " << tphase << endl;



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



/***********************************************************************************

    Direct LLL reduction

**************************************************************************************/

template<class ZT, class FT, class MatrixFT>  int
FPTuple<ZT, FT, MatrixFT>::lll(ZZ_mat<mpz_t>& C,  long alpha,
                               int lllmethod,  int sizemethod, double delta) {



  int i, j;

  int m = 1;

  ZZ_mat<ZT> L;
  L.resize(1, d);

  FP_NR<mpfr_t> t;

  for (int j = 0; j < d; j++) {
    t.mul_2si( fpv[j], alpha);
    L(0, j).set_f(t);
  }


  int found = 0;

  ZZ_mat<ZT> A;
  A.resize(m + d, d);

  // **** m=1 for the moment
  for (j = 0; j < d; j++)
    A(0, j) = L(0, j);

  for (i = 0; i < d; i++)
    A(m + i, i) = 1;


  if (lllmethod == HLLL) {



    Lattice<ZT, FT, matrix<Z_NR<ZT> >,  MatrixFT > B(A);

    B.hlll(delta);


    A = B.getbase();

  }


  else if (lllmethod == FPLLL) {

    ZZ_mat<ZT> T;
    T.resize(d, d + m);

    transpose(T, A);

    // Put the wrapper also if large examples
    lll_reduction(T, delta, 0.51, LM_FAST, FT_DEFAULT, 0);

    transpose(A, T);

  } // end FPLLL


  int foundcol = 0;

  ZZ_mat<ZT> U;

  U.resize(d, d);

  for (i = 0; i < d; i++)
    for (j = 0; j < d; j++)
      U(i, j) = A(i + 1, j);


  matprod_in(L, U);


  Z_NR<ZT> tz, tt, maxcol;

  tz.abs(L(0, 0));

  for (i = 1; i < d; i++) {

    tt.abs(L(0, i));

    if (tt.cmp(tz) == -1) {
      tz = tt;
      foundcol = i;
    }

  }


  C.resize(d, 1);
  for (j = 0; j < d; j++)
    C(j, 0) = A(m + j, foundcol);



  // cout << endl;
  // HPLLL_INFO(def + bitsize, " bits used");
  // HPLLL_INFO("Candidate relation found with bit confidence: ", -gap.exponent());

  // cout << endl << "Time lll: " << tlll << endl;
  // cout << "Time products: " << tprod << endl << endl;
  // return 1;


  // TODO: get status of LLL and assign found accordingly
  return found; // 0 here


};




} // end namespace hplll


#endif











