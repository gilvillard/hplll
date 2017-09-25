/*
Ven  3 jui 2016 15:04:49 CEST
Copyright (C) 2016      Gilles Villard

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


#include "hplll.h"

#include "wrappers.h"

#include <NTL/LLL.h>

using namespace NTL;

/* ***********************************************

          MAIN

   ********************************************** */

using namespace hplll;

int main(int argc, char *argv[])  {

  char results[] = "benchmarks_results/fhn_ntrul.results";  // ******** SPECIALIZE

  filebuf fb;
  iostream os(&fb);
  fb.open (results, ios::out);

  ZZ_mat<mpz_t> A; // For hpLLL
  ZZ_mat<mpz_t> AT, tmpmat; // fpLLL

  // ---------------------------------------------------------------------

  int k, K;

  vector<char*> s(100);

  k = 0;

  //------------


  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_80.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_140.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_200.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_240.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_280.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_320.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_360.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_400.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_440.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_480.txt", 599);
  // k += 1;


  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_520.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_560.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_600.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_640.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_680.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_720.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_800.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_900.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_1000.txt", 599);
  // k += 1;

  // s[k] = (char *) malloc(600);
  // strncpy(s[k], "collection/ntrul/ntrul_1200.txt", 599);
  // k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_181.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_227.txt", 599);
  k += 1;

s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_317.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_421.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_523.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_613.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_709.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_821.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_907.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_1061.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_1123.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_1213.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_1301.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_1427.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntrul/ntrul_1523.txt", 599);
  k += 1;

  


  //-------------

  K = k;

  double delta = 0.99;

  int run = 0;

  Timer time;

  int status = 0;

  os << endl << "FPLLL, HPLLL, NTL running times / ntrul lattice bases " << endl;   // ******** SPECIALIZE
  os << endl << "FPLLL and HPLLL hand wrapper " << endl;
  //os << endl <<  "NTL XD infinite loop for 224, with G_LLL" << endl;
  //os << endl <<  "   then NTL limited to less than 420" << endl;

  os <<         "----------------------------------------------------------------" << endl << endl;

  for (int k = 0; k < K; k++) {


    /*****************************************************************************/
    /*   i-th run  */
    /*****************************************************************************/

    run += 1;

    // Input basis
    fb.close();
    fb.open (s[k], ios::in);
    os >>  AT ;
    fb.close();
    fb.open (results, ios::app);

    int d = AT.get_rows();
    int n = AT.get_cols();

    A.resize(n, d);
    tmpmat.resize(n, d);

    transpose(A, AT);


    ZZ_mat<long> Along, ATlong;

    Along.resize(n, d);

    // -------------------

    cout << n <<  endl;

    cout << "--------------  HLLL" << endl << endl;

    {

      os << endl << "------------------------------------------------ " << endl ;

      time.start();

      matrix_cast(ATlong, AT);

      cout << endl << "HLLL calling FPLLL" << endl << endl;

      status = lll_reduction(ATlong, delta, 0.501, LM_FAST, FT_DOUBLE, 0, LLL_VERBOSE);


      transpose(Along, ATlong);

      verboseDepth = 1;

      if (status != 0) {

        cout << endl << "HLLL long + double after FPLLL" << endl << endl;

        Lattice<long, double, matrix<Z_NR<long> >,  matrix<FP_NR<double> > > B(Along, NO_TRANSFORM, DEF_REDUCTION); //* name

        status = B.hlll(delta);

        set(Along, B.getbase());


      }

      if (status != 0) {

        cout << endl << "HLLL long + long double" << endl << endl;

        Lattice<long, long double, matrix<Z_NR<long> >,  matrix<FP_NR<long double> > > B(Along, NO_TRANSFORM, DEF_REDUCTION); //* name

        status = B.hlll(delta);

        set(Along, B.getbase());


      }
      // if (n <= DIM_PREC_1) status = B.hlll(delta); //* name

      // else {
      //   hlll<long>(tmpmat, A, 0.99, false, false); // The checks are done below
      //   status = 0;
      // }

      // ***** Attemp directly without the wrapper, for the moment the heuristic stop may not be 100% sure



      if (status != 0) {

        cout << endl << "HLLL SEYSEN long + double" << endl << endl;

        Lattice<long, double, matrix<Z_NR<long> >,  matrix<FP_NR<double> > > BB(Along, NO_TRANSFORM, SEYSEN_REDUCTION);

        status = BB.hlll(delta);

        set(Along, BB.getbase());

      }

      if (status != 0) {

        cout << endl << "HLLL SEYSEN long + long double" << endl << endl;

        Lattice<long, long double, matrix<Z_NR<long> >,  matrix<FP_NR<long double> > > BBB(Along, NO_TRANSFORM, SEYSEN_REDUCTION);

        status = BBB.hlll(delta);

        set(Along, BBB.getbase());

      }

      if (status != 0)
        cout << endl << "HLLLP: REDUCTION DID NOT SUCCEED" << endl << endl;

      verboseDepth = 0;


      time.stop();


      os << "Run " << run << "  with dim = " << n << ",   delta = " << delta <<  endl << endl;
      os << "    hlll: " << time << endl ;
      time.print(os);
      os << endl;

      //if  (n <= DIM_PREC_1) matrix_cast(tmpmat, B.getbase());

      matrix_cast(tmpmat, Along);

      if (status == 0) {  // PB status non assigned with the wrapper, forced in entry
        //Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(tmpmat, NO_TRANSFORM, DEF_REDUCTION); //* names

        //T.isreduced(delta - 0.1); //* name

        double t, u, v, w;

        hplll::ratio<mpz_t>(tmpmat, t, u, v, w);


        cout << endl << ".. hplll log 2 Frobenius norm cond: " << t << endl;
        cout << ".. Average diagonal ratio: " << u << endl;
        cout << ".. Max diagonal ratio: " << v << endl;
        cout << ".. First vector quality: " << w << endl;
      }

      cout << endl;

      cout << "--------------  FPLLL HAND WRAPPER VERBOSE " << endl << endl;

      //if (n < 512) // To tune again
      {

        matrix_cast(ATlong, AT);

        time.start();

        status = lll_reduction(ATlong, delta, 0.501, LM_FAST, FT_DOUBLE, 0, LLL_VERBOSE);


        if (status != RED_SUCCESS)
          status = lll_reduction(ATlong, delta, 0.501, LM_FAST, FT_LONG_DOUBLE, 0, LLL_VERBOSE);

        if (status != RED_SUCCESS)
          status = lll_reduction(ATlong, delta, 0.501, LM_FAST, FT_DD, 0, LLL_VERBOSE);

        if (status != RED_SUCCESS)
          status = lll_reduction(ATlong, delta, 0.501, LM_PROVED, FT_MPFR, 212, LLL_VERBOSE);

        if (status != RED_SUCCESS)
          status = lll_reduction(ATlong, delta, 0.501, LM_PROVED, FT_MPFR, 424, LLL_VERBOSE);

        time.stop();


        os << "   fplll: " << time << endl << endl ;
        time.print(os);
        os << endl;

        if (status == RED_SUCCESS) {

          matrix_cast(AT, ATlong);
          transpose(tmpmat, AT);

          //Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(tmpmat, NO_TRANSFORM, DEF_REDUCTION); //* name
          //T2.isreduced(delta - 0.1); //* name

          double t, u, v, w;

          hplll::ratio<mpz_t>(tmpmat, t, u, v, w);


          cout << endl << ".. fplll log 2 Frobenius norm cond: " << t << endl;
          cout << ".. Average diagonal ratio: " << u << endl;
          cout << ".. Max diagonal ratio: " << v << endl;
          cout << ".. First vector quality: " << w << endl;

          cout << endl;

        }
        else
          cout << endl << "FPLLL Probably not reduced, mpfr stats not executed" << endl;

      }

      // cout << "--------------  NTL  " << endl << endl;

      // Mat<ZZ> BN;

      // // Input basis
      // fb.close();
      // fb.open ("tmp.txt", ios::out);
      // os <<  transpose(A) ;
      // fb.close();
      // fb.open ("tmp.txt", ios::in);
      // os >> BN;
      // fb.close();
      // system("rm tmp.txt");
      // fb.open (results, ios::app);


      // time.start();

      // if (n <= DIM_PREC_1)
      //   LLL_FP(BN, 0.99, 0, 0, 1);
      // else
      //   LLL_XD(BN, 0.99, 0, 0, 1);

      // time.stop();

      // fb.close();
      // fb.open ("tmp.txt", ios::out);
      // os <<  BN ;
      // fb.close();
      // fb.open ("tmp.txt", ios::in);
      // os >> AT;
      // fb.close();
      // system("rm tmp.txt");
      // fb.open (results, ios::app);



      // transpose(tmpmat, AT);

      // Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(tmpmat, NO_TRANSFORM, DEF_REDUCTION);
      // verboseDepth = 0;
      // T.isreduced(delta - 0.1);


      // os << "   ntl: " << time << endl << endl ;
      // time.print(os);
      // os << endl;

    }

  }// End on runs, k loop


// END
  fb.close();



  return 0;
}
