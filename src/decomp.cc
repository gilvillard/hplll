/* FGAS decomposition

Created Mar 18 jan 2011 18:10:25 CET
Main update Jeu  7 mar 2013 14:14:54 CET
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


#ifndef HPLLL_DECOMP_CC
#define HPLLL_DECOMP_CC


namespace hplll {

// ********   ATTENTION   **********************
//
// Nullity test ? Especially with longdouble ?
//
// Should be full rank on the left
//
// Attention !!! Change of precision, put back the previous one if
// needed since this is done through  mpfr_set_default_prec(setprec);
//
// Todo : proprifications, which ones and R update ?
//
//=========================================================================
//
//  Global variable dec used for switching between HJLS and PSLQ
//   dec is the dimension e.g. dec=1 one vector of the space
//   for simultaneous relations detection
//
//  Assume that the first dec columns of the fgas are, d x dec
//  are for input columns for searching relations
//  The rest of the fgas is the identity (or something else) for HJLS
//
//
// FGAS Decomposition following the paper
//
// Here for a mpfr matrix
// Precision implicitly given from the construction hence the input FGAS
//
// Should be extended for stopping at a given solution size
// Target dim : dimension of the lattice part
// Hence in the nullspace case should be lower than the input rank
//
// The transformation applied to F is U: we compute W= Transpose (U^(-1))
//              W ==> V, the name in output
//
//==========================================================================

// Difference with decompz.cc : computations ans tests on the residue Y

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> int
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::decomp(double gamma, long int targetdim, long int targetsize) {

  // Verbose
  //cout << "Direct mpfr (floating point) FGAS decomposition ...."  << endl;

  int i, j;

  vector<FP_NR<FT> >  g(d);
  g[0] = gamma;
  for (i = 1; i < d; i++) g[i].mul(g[i - 1], g[0]);

  FP_NR<FT> tmp1, tmp2;
  FP_NR<RT> r1, r2;

  int maxh = dec;

  FP_NR<FT> diagold; // For the epsilon-test
  diagold = 0.0;
  FP_NR<RT> maxYold; // For the confidence gap test
  FP_NR<RT> maxY;

  // Initial size reduction and implicit Householder
  // ***********************************************

  for (j = 0; j < n; j++) hsizereduce(j); // Implicitly computes Householder

  // Initialization for the confidence gap value
  // minimum of the diagonal R or in Y
  // *******************************************


  if (dec > 0) {

    /*maxYold.abs(Y.get(d-ldim-1,0));
      for (j=0; j<dec; j++) {
      r1.abs(Y.get(d-ldim-1,j));
      if (r1.cmp(maxYold) > 0) maxYold=r1;
      }*/

    maxYold.abs(Y.get(0, 0));
    for (i = 0; i < d - ldim; i++) {
      for (j = 0; j < dec; j++) {
        r1.abs(Y.get(i, j));
        if (r1.cmp(maxYold) > 0) maxYold = r1;
      }
    }

  }
  else {
    diagold.abs(R.get(dec, dec));
    for (i = dec + 1; i < d; i++) {
      tmp1.abs(R.get(i, i));
      if (tmp1 < diagold)
        diagold = tmp1;
    }
  }

  // Initial residue test, especially for input cases
  // with very different magnitudes, several steps
  // needed for preparation
  // -------------------------

  if (dec > 0) {
    r1 = 0.0;
    for (i = 0; i < d; i++)
      for (j = 0; j < dec; j++) {
        r2.abs(Y.get(i, j));
        if (r2.cmp(r1) > 0) r1 = r2;
      }

    // keep r1 here
    // Heuristic detection of exhaustion of precision

    if (epsilon.cmp(r1) > 0) {


      cout << " **** #tests = " << nblov << "  Anomaly: the RT precision seems exhausted" << endl;
      cout << "      Max residue component < ";
      r1.print(); cout << " for epsilon = "; epsilon.print(); cout << endl;

      return -1;
    }
  } // end heuristic on residue with dec


  // **********************************************************
  // Main loop
  // The update of R is done at the end of the work in the loop
  // **********************************************************


  while (true) {


    nblov++;

    if (((nblov % 800000) == 0) && (nblov > 0))   cout << nblov << " tests" << endl;

    // For the epsilon test Y or R
    // ***************************

    if (dec > 0) {
      maxY.abs(Y.get(d - ldim - 1, 0));
      for (j = 0; j < dec; j++) {
        r1.abs(Y.get(d - ldim - 1, j));
        if (r1.cmp(maxY) > 0) maxY = r1;
      }


      //if (maxY.cmp(epsilon) < 0) r1=0.0; // Mer  9 oct 2013 09:35:24 CEST when all relations are found
      // at the same iteration
      //else r1.div(maxY,maxYold);

      r1.div(maxY, maxYold);
      tmp1 = 1.0; // For the dummy test below

    }
    else {

      tmp1.div(R.get(d - ldim - 1, d - ldim - 1), diagold);
      tmp1.abs(tmp1);


      r1 = 1.0; // For the dummy test below
    }
    // tmp1 and r1 kept for the nextcoming test


    // New dimension test for the lattice component
    // ********************************************

    // Heuristic test,  purely input data dependent (no link with the FT précision?)
    // Anyway, problem here detected later: the "relation" will not be a relation


    if ( (r1.cmp(confidence_gap) < 0) || (tmp1.cmp(confidence_gap) < 0)) {


      // We put the non-zero column at the end for respecting the decomp structure
      if (maxh < n - ldim - 1) {

        F.colswap(maxh, n - ldim - 1);
        W.colswap(maxh, n - ldim - 1);

        if ((transf) || (dec == 0))
          U.colswap(maxh, n - ldim - 1);

        if (dec > 0)
          Y.rowswap(n - ldim - 1 - dec, maxh - dec);

        // No need of swapping in R that will be updated later
        // Note that the swap here is in the right rectangular part
      }

      ldim += 1;

    }


    // Global termination test and stopping heuristics
    // ***********************************************

    // Can be put in the while condition
    if (ldim >= targetdim) {  // Termination could be anticipated a few by testing zero results

      return (0);
    }

    // For the epsilon test, minimum of Y or check in R
    // ------------------------------------------------
    if (dec > 0) {
      maxYold = maxY;
    }
    else {
      // In case dec = 0, dec could be put to 0
      diagold.abs(R.get(dec, dec)); // Test only R.get(d-ldim-1,d-ldim-1) for efficiency?
      for (i = dec + 1; i < d - ldim - 1; i++) {
        tmp1.abs(R.get(i, i));
        if (tmp1 < diagold)
          diagold = tmp1;
      }
    }

    // Current lower bound on the relation bit size
    // k-th minima
    // --------------------------------------------
    // Redo what has been done in the previous loop for diagold
    tmp1.abs(R.get(dec, dec));
    for (j = dec + 1; j < d - ldim; j++) {
      tmp2.abs(R.get(j, j));
      if (tmp2.cmp(tmp1) > 0) tmp1 = tmp2;
    }
    tmp2 = 1.0;
    tmp1.div(tmp2, tmp1); // absolute bound > 1/max{h_ii}


    // Bit size
    Z_NR<mpz_t> res;
    set_f(res, tmp1);
    relation_bound = size_in_bits(res);

    // Stopping based on the relation bound
    // ------------------------------------

    if (relation_bound > targetsize) {
      cout << " **** No relation less than "
           << targetsize << " bits found" << endl;

      return -1;
    }

    // Current residue magnitude
    // -------------------------

    if (dec > 0) {
      r1 = 0.0;
      for (i = 0; i < d; i++)
        for (j = 0; j < dec; j++) {
          r2.abs(Y.get(i, j));
          if (r2.cmp(r1) > 0) r1 = r2;
        }



      // keep r1 here
      // Heuristic detection of exhaustion of precision
      // ----------------------------------------------


      if (epsilon.cmp(r1) > 0) {


        cout << " **** #tests = " << nblov << "  Anomaly: the RT precision seems exhausted" << endl;
        cout << "      Max residue component < ";
        r1.print(); cout << " for epsilon = "; epsilon.print(); cout << endl;

        return -1;
      }
    } // end heuristic on residue with dec
    // Heuristic on the size of U
    if ((transf) || (dec == 0)) {
      int  l = 0;
      for (i = 0; i < 2 ; i++) // 2 Instead of rowU, put hard for the heuristic check here only
        for (j = 0; j < n; j++)
          l = max(l, size_in_bits(U.get(i, j)));


      if (((unsigned int) l) > getprec_fgas() - shiftU) {

        cout << " **** #tests = " << nblov << "  Anomaly: the RT precision may be exhausted" << endl;
        cout << "      Bit size of U > " << l << endl;

        return -1;

      }

    } // end heuristic check dec==0

    // Swap position search in the left square part
    // ********************************************

    maxh = dec;
    tmp1.mul(g[0], R.get(dec, dec));
    tmp1.abs(tmp1);

    for (j = dec + 1; j < d - ldim; j++) {

      tmp2.mul(g[j], R.get(j, j));
      tmp2.abs(tmp2);

      if (tmp2 > tmp1) {
        maxh = j;  // (si plusieurs identiques ?)
        tmp1 = tmp2;
      }
    }


    // Swap ok in the left part
    // ************************

    if (maxh < d - ldim - 1) {



      F.colswap(maxh, maxh + 1);
      W.colswap(maxh, maxh + 1);


      if ((transf) || (dec == 0))
        U.colswap(maxh, maxh + 1);

      if (dec > 0)
        Y.rowswap(maxh - dec, maxh + 1 - dec);


      // Restricted to the left part, right part below in the else
      for (j = maxh; j < d - ldim; j++) hsizereduce(j); // Householder implicitly computed

    }

    // Swap position search in the right rectangular part
    // **************************************************

    else {

      // On ne fait Householder et la proprification que si on
      // rentre dans la zone de recherche

      for (j = d - ldim; j < n - ldim; j++) hsizereduce(j); // Householder implicitly computed

      maxh = d - ldim;
      tmp1.abs(R.get(d - ldim - 1, d - ldim));

      for (j = d - ldim + 1; j < n - ldim; j++) {

        tmp2.abs(R.get(d - 1 - ldim, j));

        if (tmp2 > tmp1) {
          maxh = j;
          tmp1 = tmp2;
        }

      }

      F.colswap(d - ldim - 1, maxh);
      W.colswap(d - ldim - 1, maxh);

      if ((transf) || (dec == 0))
        U.colswap(d - ldim - 1, maxh);

      if (dec > 0)
        Y.rowswap(d - ldim - 1 - dec, maxh - dec);

      // Needed: this column can be used for swap if non-zero diag at next phase
      // Not reduced everywhere for not reducing w.r.t. a zero (very small) entry
      //  the rest will be reduced at next phase in the right part if non-zero diag

      hsizereduce(d - ldim - 1);

    }  // End else right part

  } // End main iteration loop

  return 0;

}


/* -------------------------------------------------------------------------
   Size reduction

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm
   ------------------------------------------------------------------------- */
// exactly same as in decompz.cc but FP_NR<RT> xr; necessary here
// the global shift dec variable is used

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> int
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::hsizereduce(int kappa) {


  FP_NR<FT> approx;

  approx = 0.01;

  FP_NR<FT> x, t, tmpfp;
  Z_NR<ZT>  xz, tmpz;

  FP_NR<RT> xr;

  int i, w = 0;

  bool nonstop = 1;
  bool somedone = 0;

  int elim_index;


  householder(kappa);

  while (nonstop) {

    w++;

    somedone = 0;

    elim_index = min(kappa - 1, d - 1 - ldim);


    for (i = elim_index; i > -1; i--) {

      x.div(R.get(i, kappa), R.get(i, i));

      x.rnd(x);

      if (x.sgn() != 0) {

        set_f(xz, x);

        if (xz == -1) {
          somedone = 1;

          R.addcol(kappa, i, i);

          F.addcol(kappa, i, d);
          W.subcol(i, kappa, n);

          if ((transf) || (dec == 0))
            U.addcol(kappa, i, rowU);

          if ( (i >= dec) && (dec > 0))
            Y.subrow(i - dec, kappa - dec);

        }
        /* ----------------------------------------------- */
        else if (xz == 1) {
          somedone = 1;

          R.subcol(kappa, i, i);

          F.subcol(kappa, i, d);
          W.addcol(i, kappa, n);

          if ((transf) || (dec == 0))
            U.subcol(kappa, i, rowU);

          if ( (i >= dec) && (dec > 0))
            Y.addrow(i - dec, kappa - dec);

        }
        /* ----------------------------------------------- */
        else {
          somedone = 1;

          set_z(xr, xz);

          R.submulcol(kappa, i, x, i);
          F.submulcol(kappa, i, xr, d);

          /*if (x.exponent() < 31) {

            int ii;
            ii=xz.get_si();

            if (ii > 0)
              vector_addmul_ui(W.getcol(i),W.getcol(kappa),ii,n);
            else
              vector_submul_ui(W.getcol(i),W.getcol(kappa),-ii,n);
          }

          else {*/

          W.addmulcol(i, kappa, xz, n);
          //}

          if ((transf) || (dec == 0))
            U.submulcol(kappa, i, xz, rowU);

          if ( (i >= dec) && (dec > 0))
            Y.addmulrow(i - dec, kappa - dec, xr);


        }
        /* ----------------------------------------------- */

      } // End non zero combination

    } // End loop through the column


    if (somedone) {

      t.mul(approx, normB2[kappa]);

      householder(kappa);  // Implicit update of normB2

      nonstop = (normB2[kappa] < t);  // Some decrease ?

    }
    else {
      nonstop = 0;
    }
  } // end while


  return 0;

};



/* --------------------------------------------- */
/* Householder on column kappa                   */
/*     partially : d-ldim x  n                   */
/* --------------------------------------------- */
// exactly same as in decompz.cc

// Change w.r.t. hlll _up and _v and we consider ldim


template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> int
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::householder(int kappa) {

  int k;
  FP_NR<FT> nrtmp, s, w;


  // Left non-singular part
  // ----------------------
  if (kappa < d - ldim)  {

    R.setcol(kappa, F.getcol(kappa), 0, d);

    fp_norm_sq(normB2[kappa], R.getcol(kappa), d);

    for (k = 0; k < kappa; k++) {

      scalarprod(nrtmp, V.getcol(k, k), R.getcol(kappa, k), d - k);

      R.fmasub(kappa, k, R.getcol(kappa, k), V.getcol(k, k), nrtmp, d - k);

    }

    w = R.get(kappa, kappa);

    if (w >= 0) {
      fp_norm(s, R.getcol(kappa, kappa), d - kappa);
      nrtmp.neg(s);
      R.set(kappa, kappa, nrtmp);
    }
    else {
      fp_norm(nrtmp, R.getcol(kappa, kappa), d - kappa);
      R.set(kappa, kappa, nrtmp);
      s.neg(nrtmp);
    }

    w.add(w, s);
    s.mul(s, w);
    s.sqrt(s);

    V.div(kappa, kappa + 1, R.getcol(kappa, kappa + 1), s, d - kappa - 1);
    nrtmp.div(w, s);
    V.set(kappa, kappa, nrtmp);

    for (int i = kappa + 1; i < d; i++) R.set(i, kappa, 0.0);

  }  // end of kappa in left part

  else {

    // Non square part
    // ---------------


    R.setcol(kappa, F.getcol(kappa), 0, d);

    fp_norm_sq(normB2[kappa], R.getcol(kappa), d);

    for (k = 0; k < d - ldim; k++) { // Introduction of ldim

      scalarprod(nrtmp, V.getcol(k, k), R.getcol(kappa, k), d - k);

      R.fmasub(kappa, k, R.getcol(kappa, k), V.getcol(k, k), nrtmp, d - k);

    }
  } // end of kappa in right part


  return 0;
}


/* --------------------------------------------- */
/* Complete Householder  d x  n                  */
/* --------------------------------------------- */
// exactly same as in decompz.cc

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> int
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::householder() {

  int k, kappa;
  FP_NR<FT> nrtmp, s, w;

  // Full rank left square part
  // --------------------------

  for (kappa = 0; kappa < d; kappa++) {

    R.setcol(kappa, F.getcol(kappa), 0, d);

    for (k = 0; k < kappa; k++) {

      scalarprod(nrtmp, V.getcol(k, k), R.getcol(kappa, k), d - k);
      R.fmasub(kappa, k, R.getcol(kappa, k), V.getcol(k, k), nrtmp, d - k);

    }

    w = R.get(kappa, kappa);

    if (w >= 0) {
      fp_norm(s, R.getcol(kappa, kappa), d - kappa);
      nrtmp.neg(s);
      R.set(kappa, kappa, nrtmp);
    }
    else {
      fp_norm(nrtmp, R.getcol(kappa, kappa), d - kappa);
      R.set(kappa, kappa, nrtmp);

      s.neg(nrtmp);
    }

    w.add(w, s);
    s.mul(s, w);
    s.sqrt(s);

    V.div(kappa, kappa + 1, R.getcol(kappa, kappa + 1), s, d - kappa - 1);
    nrtmp.div(w, s);
    V.set(kappa, kappa, nrtmp);

    for (int i = kappa + 1; i < d; i++) R.set(i, kappa, 0.0);

  }  // end left part

  // Non square part
  // ---------------

  for (kappa = d; kappa < n; kappa++) {

    R.setcol(kappa, F.getcol(kappa), 0, d);

    for (k = 0; k < d; k++) {

      scalarprod(nrtmp, V.getcol(k, k), R.getcol(kappa, k), d - k);
      R.fmasub(kappa, k, R.getcol(kappa, k), V.getcol(k, k), nrtmp, d - k);

    }
  } // end right part


  return 0;
}


// ********************************
// VARIOUS   ACCESS CONSTRUCTION
// ********************************

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
inline unsigned int Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::getprec_fgas() {

  // Actually the global ambiant mpfr prec (see also setprec_internal)
  // FP_NR getprec is global, not the one of the element with mpfr
  return (F.get(0, 0)).getprec();

}



// For the zero test or exhausted RT precision / to tune
// Modify the global fgas variable epsilon
// -----------------------------------------------------

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
inline unsigned int Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::shift_epsilon(int bits) {

  int oldprec = getprec_fgas();

  shiftepsilon = bits;

  mpfr_set_default_prec(100); // !!!!!!!
  FP_NR<mpfr_t> r1, r2, ee;
  double prec;

  // !! bits and oldprec must be int
  prec = (double) (shiftepsilon - oldprec); // heuritic to check

  r1 = prec;
  r2 = 2.0;
  r2.log(r2);
  r1.mul(r1, r2);
  ee.exponential(r1);
  set_mpfr(epsilon, ee);
  mpfr_set_default_prec(oldprec);


  return (oldprec);

}

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
inline unsigned int Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::shift_testU(int bits) {



  shiftU = bits;

  return (shiftU);

}



template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
inline  int Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::set_confidence_gap(double ratio) {

  confidence_gap = ratio;

  return 1;
}


template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
inline ZZ_mat<ZT> Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::getU() {

  ZZ_mat<ZT> UU(n, n);


  if (transf) {

    UU.resize(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UU.Set(i, j, U(i, j)); // reprendre boucle sur les colonnes

    return UU;
  }
  else {

    cout << "*** Access error: U has not been computed" << endl;
    return UU;

  }
}


template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
inline ZZ_mat<ZT> Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::getV() {

  ZZ_mat<ZT> VV(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) VV.Set(i, j, W(i, j)); // reprendre boucle sur les colonnes

  return VV;
}



// Setprec for the internal FT elements, not for RT outside
// --------------------------------------------------------

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> unsigned int
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::setprec_internal(long prec) {

  unsigned oldglobal;
  oldglobal = getprec_fgas(); // should be put befor the mpfr change (see getprec_fgas)

  unsigned oldprec;
  oldprec = (R.get(0, 0)).getprec();

  mpfr_set_default_prec(prec);
  R.clear();
  R.resize(d, n);
  unsigned newprec;
  newprec = (R.get(0, 0)).getprec();
  if (newprec == oldprec) cout << "Warning: in function setprec_internal decomp, the change of precision has no effect" << endl;

  V.clear();
  V.resize(d, d);

  normB2.clear();
  normB2.resize(n);

  // TODO change epsilon !!!

  mpfr_set_default_prec(oldglobal);
  return oldprec;

}


template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
matrix<FP_NR<FT> > Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::getR()
{
  matrix<FP_NR<FT> >  RR(d, n);
  FP_NR<FT> tmp;

  for (int i = 0; i < d; i++)
    for (int j = i; j < n; j++) {
      tmp = R.get(i, j); // cf l'absence de const dans nr.cpp Set / Exp
      RR.set(i, j, tmp); // reprendre boucle sur les colonnes

    }
  for (int i = 0; i < d; i++)
    for (int j = 0; j < i; j++) RR.set(i, j, 0.0);

  return RR;
}


template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
MatrixRT Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::getfgas()
{
  return F;
}





// Initialization and construction
// -------------------------------

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT> void
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::init(int d, int n, bool transform, long setprec, long int inputdec) {

  dec = inputdec;

  tmpcompt = 0;

  nblov = 0;

  relation_bound = 0;

  ldim = 0;

  // Inverse transpose of the transformation matrix (hence possibly nullspace)
  W.resize(n, n);
  for (int i = 0; i < n; i++) W(i, i) = 1;

  // !!! Put back the old precision if needed outside
  // Pas forcément MPFR !!!
  mpfr_set_default_prec(setprec);

  // Verbose
  //cout << endl << "Decomp precision for mpfr: " <<  mpfr_get_default_prec() << " bits" << endl ;

  F.resize(d, n);

  transf = transform;
  if (transf)
    rowU = n;
  else
    rowU = 2;

  // Transformation U and for the heuritic checks
  if ((transf) || (dec == 0)) {
    rowU = n;
    U.resize(rowU, n);
    for (int i = 0; i < rowU; i++) U(i, i) = 1;
  }

  if (dec > 0)
    Y.resize(d, dec);


  R.resize(d, n);

  V.resize(d, d);

  normB2.resize(n);

  // For the zero test or exhausted RT precision / to tune
  // -----------------------------------------------------


  shiftU = 0;

  shiftepsilon = n + 200; // heuritic to check

  mpfr_set_default_prec(100); // Temporary
  FP_NR<mpfr_t> r1, r2, ee;
  double prec;

  prec = (double) shiftepsilon - setprec;

  r1 = prec;
  r2 = 2.0;
  r2.log(r2);
  r1.mul(r1, r2);
  ee.exponential(r1);

  set_mpfr(epsilon, ee);
  mpfr_set_default_prec(setprec);

  // Ratio non zero / zero
  confidence_gap = 1e-9;


}


// Construction
template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::Fgas(MatrixRT Finput, bool transform, long setprec, long int inputdec) {

  d = Finput.getRows();
  n = Finput.getCols();

  //n=d-m;

  init(d, n, transform, setprec, inputdec);

  int i, j;

  for (i = 0; i < d; i++)
    for (j = 0; j < n; j++) {

      F.set(i, j, Finput.get(i, j));
    }

  if (dec > 0) {
    for (i = 0; i < d; i++)
      for (j = 0; j < dec; j++)
        Y.set(i, j, Finput.get(i, j));

  }

}

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::Fgas(MatrixRT Finput, MatrixRT res, bool transform, long setprec, long int inputdec) {

  d = Finput.getRows();
  n = Finput.getCols();

  //n=d-m;

  init(d, n, transform, setprec, inputdec);

  int i, j;

  for (i = 0; i < d; i++)
    for (j = 0; j < n; j++) {
      F.set(i, j, Finput.get(i, j));
    }

  if (dec > 0) {
    for (i = 0; i < d; i++)
      for (j = 0; j < dec; j++)
        Y.set(i, j, res.get(i, j));
  }

}



template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>  void
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::set(MatrixRT Finput) {

  int i, j;

  for (i = 0; i < d; i++)
    for (j = 0; j < n; j++) {
      F.set(i, j, Finput.get(i, j));
    }


  if (dec > 0) {
    for (i = 0; i < d; i++)
      for (j = 0; j < dec; j++)
        Y.set(i, j, Finput.get(i, j));
  }

  if ((transf) || (dec == 0)) {
    U.resize(rowU, n);
    for (int i = 0; i < rowU; i++) U(i, i) = 1;
  }

  W.resize(n, n);
  for (int i = 0; i < n; i++) W(i, i) = 1;

}

template<class RT, class ZT, class FT, class MatrixRT, class MatrixZT, class MatrixFT>  void
Fgas<RT, ZT, FT, MatrixRT, MatrixZT, MatrixFT>::set(MatrixRT Finput, MatrixRT res) {

  int i, j;

  for (i = 0; i < d; i++)
    for (j = 0; j < n; j++) {
      F.set(i, j, Finput.get(i, j));
    }

  if (dec > 0) {
    for (i = 0; i < d; i++)
      for (j = 0; j < dec; j++)
        Y.set(i, j, res.get(i, j));
  }

  if ((transf) || (dec == 0)) {
    U.resize(rowU, n);
    for (int i = 0; i < rowU; i++) U(i, i) = 1;
  }

  W.resize(n, n);
  for (int i = 0; i < n; i++) W(i, i) = 1;

}

} // end namespace hplll


#endif


