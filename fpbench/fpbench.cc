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


// ***********************************************


int main(int argc, char *argv[])  {


	// typedef double FT;


	// cout << endl << endl <<  "               FP_NR < double >          " << endl;


	// bench<FT>(vaxpy_in<FT>, 800);
	// bench<FT>(dotproduct<FT>, 800);
	// bench<FT>(vadd<FT>, 800);
	// bench<FT>(vdiv<FT>, 800);


	// typedef long double FT;


	// cout << endl << endl <<  "               FP_NR < long double >          " << endl;


	// bench<FT>(vaxpy_in<FT>, 400);
	// bench<FT>(dotproduct<FT>, 400);
	// bench<FT>(vadd<FT>, 400);
	// bench<FT>(vdiv<FT>, 400);



	cout << endl << endl <<  "               FP_NR < dd_real >          " << endl;


	typedef dd_real FT;

	bench<FT>(vaxpy_in<FT>, 400);
	bench<FT>(dotproduct<FT>, 400);
	bench<FT>(vadd<FT>, 400);
	bench<FT>(vdiv<FT>, 400);


	return 0;
}









