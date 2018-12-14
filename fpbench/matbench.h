/*

Created Mar 11 déc 2018 17:27:16 CET
Copyright (C) 2018      Gilles Villard

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

#ifndef HPLLL_MATBENCH_H
#define HPLLL_MATBENCH_H

#include "hplll.h"

namespace hplll {

// ***********************************************


// Square matrix
template<class FT> void mmul(string& tag, double& time, ZZ_mat<mpz_t> A) {


	clock_t fl_start, fl_end;

	int n = A.get_rows();

	ZZ_mat<mpz_t> B;
	B.resize(n, n);
	set(B, A);

	ZZ_mat<mpz_t> C;
	C.resize(n, n);

	fl_start = clock();

	matprod(C, B, A);

	fl_end = clock();

	time = (difftime(fl_end, fl_start) / CLOCKS_PER_SEC);

	tag = "matrix product";


}


// Square matrix
template<class FT> void householder(string& tag, double& time, ZZ_mat<mpz_t> A0) {


	clock_t fl_start, fl_end;

	Lattice<mpz_t, FT, matrix<Z_NR<mpz_t> >, matrix<FP_NR<FT> > > A(A0);

	fl_start = clock();

	A.householder();

	fl_end = clock();

	time = (difftime(fl_end, fl_start) / CLOCKS_PER_SEC);

	tag = "householder";


}

template<class FT> void gso(string& tag, double& time, ZZ_mat<mpz_t> A) {


	clock_t fl_start, fl_end;

	transpose(A, A);

	ZZ_mat<mpz_t> uz;
	ZZ_mat<mpz_t> u_invZ;

	int gso_flags = 0;
	gso_flags |= GSO_ROW_EXPO; 

	MatGSO<Z_NR<mpz_t>, FP_NR<FT> > m_gso(A, uz, u_invZ, gso_flags);

	fl_start = clock();

	m_gso.update_gso();

	fl_end = clock();

	// int n = A.get_rows();
	// print2maple(m_gso.get_r_matrix(), n, n);

	time = (difftime(fl_end, fl_start) / CLOCKS_PER_SEC);

	tag = "fplll gso";


}





// *********************************************


template<class FT> void bench(double& time, ZZ_mat<mpz_t> A, const function < void(string&, double&, ZZ_mat<mpz_t>) > &g) {

	string tag;


	int n = A.get_rows();

	g(tag, time, A);

	cout << endl << "-------------------   " << tag << "  dimension " << n << "   --------------------------" << endl;

	cout << "Time: " << time << endl;

	cout << endl;
}


} // end namespace hplll

#endif



