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

	ZZ_mat<mpz_t> C;

	ZZ_mat<mpz_t> AT;

	// ---------------------------------------------------------------------

	filebuf fb;
	iostream os(&fb);

	int n, d;
	double delta = 0.99;

	//command_line_basis(A, n, d, delta, argc, argv);

	fb.open ("basis.txt", ios::in);
	os >>  AT ;
	fb.close();

	d  = AT.get_rows();
	n = AT.get_cols();

	A.resize(n, d);

	transpose(A, AT);


	Timer th, tf;

	// HLLL ------------------------------------------
	cout << "--------------  HPLLL WRAPPER" << endl << endl;

	int status;

	ZZ_mat<long> Along;
	matrix_cast(Along, A);

	//Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A, NO_TRANSFORM, DEF_REDUCTION);
	//Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A,NO_TRANSFORM,DEF_REDUCTION);
	Lattice<long, double, matrix<Z_NR<long> >,  matrix<FP_NR<double> > >  B(Along, NO_TRANSFORM, DEF_REDUCTION);

	verboseDepth = 1;
	th.start();
	status = B.hlll(delta);
	th.stop();


	matrix_cast(A, B.getbase());

	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > TB(A, NO_TRANSFORM, DEF_REDUCTION);
	verboseDepth = 0;
	TB.isreduced(delta - 0.1);

	//th=hlll<mpz_t>(C, A, 0.99, true, false);

	//th=hlll<__int128_t>(C, A, 0.99, true,true);
	//hlll<long>(C, A, 0.99, false, true);

	time.start();

	tf.stop();

	L.relation(C, alpha, 20, 20, 40, FPLLL);
	//L.relation(C, alpha, 30, 400, -1, FPLLL);   // -1 for bits only with mpz_t
	//L.lll(C, 12220);

	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A, NO_TRANSFORM, DEF_REDUCTION);
	verboseDepth = 0;
	T.isreduced(delta - 0.1);

	//  cout << "-----------------------" << endl;

	cout << endl;

	cout << "HLLL: " << th << endl;

	cout << "FPLLL :" << tf << endl;


	return 0;
}
