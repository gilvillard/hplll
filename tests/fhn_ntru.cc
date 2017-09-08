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

  char results[] = "benchmarks_results/fhn_ntru.results";  // ******** SPECIALIZE

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


  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntru/ntru_80_53.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntru/ntru_120_80.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntru/ntru_160_107.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntru/ntru_180_120.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntru/ntru_200_133.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntru/ntru_220_147.txt", 599);
  k += 1;

  s[k] = (char *) malloc(600);
  strncpy(s[k], "collection/ntru/ntru_240_160.txt", 599);
  k += 1;

  cout << s[0] << endl;

  //-------------

  K = k;

  double delta = 0.99;

  int run = 0;

  Timer time;

  int status = 0;

  os << endl << "FPLLL, HPLLL, NTL running times / ntru lattice bases " << endl;   // ******** SPECIALIZE
  os << endl << "HPLLL hand wrapper " << endl;
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


    cout << n << "   " << d << endl;

    cout << "--------------  HLLL" << endl << endl;

    {

      os << endl << "------------------------------------------------ " << endl ;

      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  B(A, NO_TRANSFORM, DEF_REDUCTION); //* name

      time.start();

      verboseDepth = 1;

      // if (n <= DIM_PREC_1) status = B.hlll(delta); //* name

      // else hlll<mpz_t>(tmpmat, A, 0.99, true, true);


      // ***** Attemp directly without the wrapper, for the moment the heuristic stop may not be 100% sure

      status = B.hlll(delta);

      set(tmpmat, B.getbase());

      if (status != 0) {

        Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  BB(tmpmat, NO_TRANSFORM, SEYSEN_REDUCTION);

        status = BB.hlll(delta);

        set(tmpmat, BB.getbase());

      }

      if (status != 0) {

        Lattice<mpz_t, ldpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<long double, ldpe_t> >  BBB(tmpmat, NO_TRANSFORM, SEYSEN_REDUCTION);

        status = BBB.hlll(delta);

        set(tmpmat, BBB.getbase());

      }

      verboseDepth = 0;
      time.stop();


      // if (n <= DIM_PREC_1)
      //   matrix_cast(tmpmat, B.getbase());

      os << "Run " << run << "  with n,d = " << n << "  " << d << ",    delta = " << delta <<  endl << endl;
      os << "    hlll: " << time << endl ;
      time.print(os);
      os << endl;


      if (status == 0) {
        Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(tmpmat, NO_TRANSFORM, DEF_REDUCTION); //* names

        T.isreduced(delta - 0.1); //* name

        double t, u, v, w;

        hplll::ratio<mpz_t>(tmpmat, t, u, v, w);

        cout << endl << ".. hplll log 2 Frobenius norm cond: " << t << endl;
        cout << ".. Average diagonal ratio: " << u << endl;
        cout << ".. Max diagonal ratio: " << v << endl;
        cout << ".. First vector quality: " << w << endl;
      }
      cout << endl;

      cout << "--------------  FPLLL WRAPPER VERBOSE " << endl << endl;

      time.start();
      lll_reduction(AT, delta, 0.501, LM_WRAPPER, FT_DEFAULT, 0, LLL_VERBOSE);
      time.stop();


      os << "   fplll: " << time << endl << endl ;
      time.print(os);
      os << endl;

      transpose(A, AT);
      Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(A, NO_TRANSFORM, DEF_REDUCTION); //* name
      T2.isreduced(delta - 0.1); //* name

      double t, u, v, w;

      hplll::ratio<mpz_t>(A, t, u, v, w);

      cout << endl << ".. fplll log 2 Frobenius norm cond: " << t << endl;
      cout << ".. Average diagonal ratio: " << u << endl;
      cout << ".. Max diagonal ratio: " << v << endl;
      cout << ".. First vector quality: " << w << endl;

      // cout << "--------------  NTL  " << endl << endl;

      // Mat<ZZ> BN;

      // // Input basis
      // fb.close();
      // fb.open (s[k], ios::in);
      // os >>  BN ;
      // fb.close();
      // fb.open (results, ios::app);

      // time.start();

      // if (n < 224)
      //   LLL_XD(BN, 0.99, 0, 0, 1);
      // else if  (n < 420)
      //   G_LLL_XD(BN, 0.99, 0, 0, 1);

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

      // transpose(A, AT);

      // Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A, NO_TRANSFORM, DEF_REDUCTION);
      // verboseDepth = 0;
      // if  (n < 420)
      //   T.isreduced(delta - 0.1);


      // os << "   ntl: " << time << endl << endl ;
      // time.print(os);
      // os << endl;

    }

  }// End on runs, k loop


  // END
  fb.close();



  return 0;
}
