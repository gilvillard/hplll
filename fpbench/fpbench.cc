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



#include "fpbench.h"

using namespace hplll;

#include "nr_FP_float128.inl"

// ***********************************************

template<class FT> void go(vector<double>& t, int rounds) {

	bench<FT>(t[0], t[1], vaxpy_in<FT>, rounds);
	bench<FT>(t[2], t[3], dotproduct<FT>, rounds);
	bench<FT>(t[4], t[5], vadd<FT>, rounds);
	bench<FT>(t[6], t[7], vdiv<FT>, rounds);

	cout << "***************************" << endl;
	printf("Table:  %.1f    %.1f    %.1f    %.1f    %.1f    %.1f \n \n ", t[0], t[1], t[2], t[3], t[5], t[7]);

}


int main(int argc, char *argv[])  {


	vector<double> time;
	time.resize(20);


	// --------------------------

	/*typedef double FT;
	cout << endl << endl <<  "               FP_NR < double >          " << endl;

	go<FT>(time, 1200);*/


	/*typedef long double FT;
	cout << endl << endl <<  "               FP_NR < long double >          " << endl;

	go<FT>(time, 1200);*/


	/*typedef dd_real FT;
	cout << endl << endl <<  "               FP_NR < dd_real >          " << endl;

	go<FT>(time, 1200);*/


	/*mpfr_set_default_prec(106);

	typedef mpfr_t FT;
	cout << endl << endl <<  "               FP_NR < mpfr_t >  106         " << endl;

	go<FT>(time, 600);
	*/

	// mpfr_set_default_prec(212);

	// typedef mpfr_t FT;
	// cout << endl << endl <<  "               FP_NR < mpfr_t >  212        " << endl;

	// go<FT>(time, 600);


	// typedef qd_real FT;
	// cout << endl << endl <<  "               FP_NR < qd_real>         " << endl;

	// go<FT>(time, 600);


	typedef __float128 FT;
	cout << endl << endl <<  "               FP_NR < __float128 >         " << endl;

	go<FT>(time, 600);

	

	return 0;
}









