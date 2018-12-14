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



#include "matbench.h"

using namespace hplll;

#include "nr_FP_float128.inl"


// ***********************************************

template<class FT> void go(int n, int nbbits) {

	vector<double> t;
	t.resize(4);

	ZZ_mat<mpz_t> A;

	A.resize(n, n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n ; j++)
			A(i, j).randb(nbbits);


	bench<double>(t[0], A, householder<double>);

	bench<FT>(t[1], A, householder<FT>);

	bench<double>(t[2], A, gso<double>);

	bench<FT>(t[3], A, gso<FT>);

	cout << "***************************" << endl;
	printf("Table:  %.3f    %.3f    %.1f    %.3f    %.3f    %.1f \n \n ", t[0], t[1], t[1] / t[0], t[2], t[3], t[3] / t[2]);

}


int main(int argc, char *argv[])  {


	// -------------------------------------------------------------------------------

	int n = 200;

	cout << endl << endl <<  "               FP_NR < double >        " << endl;
	go<double>(n, 20);


	cout << endl << endl <<  "               FP_NR < long double >        " << endl;
	go<long double>(n, 20);


	cout << endl << endl <<  "               FP_NR < dd_real >          " << endl;
	go<dd_real>(n, 20);


	mpfr_set_default_prec(106);
	cout << endl << endl <<  "               FP_NR < mpfr_t >  106         " << endl;
	go<mpfr_t>(n, 20);


	// cout << endl << endl <<  "               FP_NR < __float128 >         " << endl;
	// go<__float128>(n,20);

	mpfr_set_default_prec(212);
	cout << endl << endl <<  "               FP_NR < mpfr_t >  212        " << endl;
	go<mpfr_t>(n, 20);


	cout << endl << endl <<  "               FP_NR < qd_real>         " << endl;
	go<qd_real>(n, 20);


	// -------------------------------------------------------------------------------

	n = 800;

	cout << endl << endl <<  "               FP_NR < double >        " << endl;
	go<double>(n, 20);


	cout << endl << endl <<  "               FP_NR < long double >        " << endl;
	go<long double>(n, 20);


	cout << endl << endl <<  "               FP_NR < dd_real >          " << endl;
	go<dd_real>(n, 20);


	mpfr_set_default_prec(106);
	cout << endl << endl <<  "               FP_NR < mpfr_t >  106         " << endl;
	go<mpfr_t>(n, 20);


	// cout << endl << endl <<  "               FP_NR < __float128 >         " << endl;
	// go<__float128>(n,20);

	mpfr_set_default_prec(212);
	cout << endl << endl <<  "               FP_NR < mpfr_t >  212        " << endl;
	go<mpfr_t>(n, 20);


	cout << endl << endl <<  "               FP_NR < qd_real>         " << endl;
	go<qd_real>(n, 20);


	return 0;
}









