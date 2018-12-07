/*

Created Dim  7 avr 2013 16:54:03 CEST
Copyright (C) 2013-2016      Gilles Villard

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

/* ***********************************************

          MAIN

   ********************************************** */

using namespace hplll;

int main(int argc, char *argv[])  {


	ZZ_mat<mpz_t> A;
	ZZ_mat<mpz_t> L;

	ZZ_mat<mpz_t> Along;

	ZZ_mat<mpz_t> AT;

	// ---------------------------------------------------------------------

	filebuf fb;
	iostream os(&fb);

	int n, d;
	double delta = 0.99;

	command_line_basis(A, n, d, delta, argc, argv);

	// fb.open ("basis.txt", ios::in);
	// os >>  AT ;
	// fb.close();

	// d  = AT.get_rows();
	// n = AT.get_cols();

	// A.resize(n, d);

	// transpose(A, AT);

	AT.resize(d, n);

	Timer tp, ts;


	//matrix_cast(Along, A);

	//SLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A, 4, NO_TRANSFORM, DEF_REDUCTION);
	//SLattice<long, double, matrix<Z_NR<long> >, matrix<FP_NR<double> > > B(Along, 4, TRANSFORM, DEF_REDUCTION);


<<<<<<< HEAD
	int S = 32;
=======
>>>>>>> b9cf7a0f0e977249d7972a607a2592cf1ef594a1

	int S = 8;

	SLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  B(A, S, NO_TRANSFORM, DEF_REDUCTION);

	tp.clear();
	tp.start();

	B.phouseholder(S);

	tp.stop();

	cout << "pLLL: " << tp << endl;

	tp.clear();
	tp.start();

	B.householder(d);

	tp.stop();


	cout << "LLL :" << tp << endl;


	// tp.start();

	// //B.hlll(delta, 2, 2);

	// slll_wrap<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > (L, A, 30, 4, delta, DEF_REDUCTION);

	// tp.stop();


	// Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A);

	// ts.start();

	// C.hlll(delta);

	// ts.stop();


	// cout << endl << "-----------------------" << endl;

	// cout << endl;

	// cout << "SLLL: " << tp << endl;

	// cout << "HPLLL :" << ts << endl;


	return 0;
}
