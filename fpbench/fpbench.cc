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


	double time[20];

	// typedef double FT;
	// cout << endl << endl <<  "               FP_NR < double >          " << endl;

	// bench<FT>(time[0], time[1], vaxpy_in<FT>, 1200);
	// bench<FT>(time[2], time[3], dotproduct<FT>, 1200);
	// bench<FT>(time[4], time[5], vadd<FT>, 1200);
	// bench<FT>(time[6], time[7], vdiv<FT>, 1200);

	// cout << "Table results: " <<  time[0] << "    " <<  time[1] << "    " <<  time[2] << "    " <<  time[3] 
	//      << "    " <<  time[5] << "    " <<  time[7] << "    " << endl;


	typedef double FT;
	cout << endl << endl <<  "               FP_NR < long double >          " << endl;

	bench<FT>(time[0], time[1], vaxpy_in<FT>, 1200);
	bench<FT>(time[2], time[3], dotproduct<FT>, 1200);
	bench<FT>(time[4], time[5], vadd<FT>, 1200);
	bench<FT>(time[6], time[7], vdiv<FT>, 1200);

	cout << "Table results: " <<  time[0] << "    " <<  time[1] << "    " <<  time[2] << "    " <<  time[3] 
	     << "    " <<  time[5] << "    " <<  time[7] << "    " << endl;



	//cout << endl << endl <<  "               FP_NR < dd_real >          " << endl;


	//typedef dd_real FT;




	return 0;
}









