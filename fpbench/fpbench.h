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

#ifndef HPLLL_FPBENCH_H
#define HPLLL_FPBENCH_H

#include "hplll.h"


namespace hplll {

// ***********************************************


template<class FT> void vadd(string& tag, int& count, int n, FP_NR<FT>& r, vector<FP_NR<FT> >& va, vector<FP_NR<FT> >& vb) {

	r = 0.0;

	for (int i = 0; i < n; i++) {

		va[i].add(va[i], vb[i]);

	}

	tag = "vector add";

	count = 1;
}


template<class FT> void dotproduct(string& tag, int& count, int n, FP_NR<FT>& r, vector<FP_NR<FT> >& va, vector<FP_NR<FT> >& vb) {

	r = 0.0;

	for (int i = 0; i < n; i++) {

		r.addmul(va[i], vb[i]);

	}

	tag = "dot product";

	count = 2;
}


template<class FT> void vaxpy_in(string& tag, int& count, int n, FP_NR<FT>& r, vector<FP_NR<FT> >& va, vector<FP_NR<FT> >& vb) {

	r = 0.0;

	FP_NR<FT> t = va[0];

	for (int i = 0; i < n; i++) {

		va[i].addmul(va[i], t);

	}

	tag = "vaxpy in";

	count = 2;
}

template<class FT> void vdiv(string& tag, int& count, int n, FP_NR<FT>& r, vector<FP_NR<FT> >& va, vector<FP_NR<FT> >& vb) {

	r = 0.0;

	for (int i = 0; i < n; i++) {

		va[i].div(va[i], vb[i]);

	}

	tag = "vdiv";

	count = 1;
}

// *********************************************


template<class FT> void bench(double& t1, double& t2, const function<void(string&, int&, int, FP_NR<FT>&, vector<FP_NR<FT> >& , vector<FP_NR<FT> >&)> &g, int nbtrials) {


	t1 = 0.0;
	t2 = 0.0;

	FP_NR<FT> f = 0.0;

	int prec = f.get_prec();

	int nsize = 2;

	vector<int> n;
	n.resize(nsize);

	n[0] = 400;
	n[1] = 10000;
	//n[2] = 100000;


	vector<FP_NR<FT> > va;
	vector<FP_NR<FT> > vb;

	FP_NR<FT> r;

	double fpr = 0.0;

	int count;

	clock_t fl_start, fl_end;


	// Measure for pure doubles
	// ------------------------


	vector<double> r1;

	int fpdim = 10000;
	int nbfp = 8000;

	r1.resize(fpdim);

	double fp_time = 0.0;

	srand(time(NULL));

	for (int K = 0; K < nbfp ; K++) {

		double c = rand();

		for (int i = 0; i < fpdim; i++)
			r1[i] = rand();

		// Chrono
		// ------

		fl_start = clock();

		for (int i = 0; i < fpdim; i++)
			r1[i] += c * r1[i];

		fl_end = clock();

		fp_time += (difftime(fl_end, fl_start) / CLOCKS_PER_SEC);

		// use res
		for (int i = 0; i < fpdim; i++)
			fpr += r1[i];

	}

	// use res
	cout << endl << fpr << endl;

	fp_time = fp_time / nbfp;

	double dflops = 2 * ((double) fpdim) / fp_time;

	cout << endl << "***** dGflops    [" << fpdim << "]: " << dflops << endl;


	// Measure for float type
	// ----------------------


	for (int l = 0; l < nsize; l++) {

		va.resize(n[l]);
		vb.resize(n[l]);


		Z_NR<mpz_t> r0;

		srand(time(NULL));
		r0.randb(prec);


		string tag;

		double fpnr_time = 0.0;


		for (int K = 0; K < nbtrials; K++) {

			for (int i = 0; i < n[l]; i++) {

				va[i] = 0.0;
				vb[i] = 0.0;

				for (int k = prec; k > 0; k /= 53) {

					va[i].mul_2si(va[i], 53);
					va[i].add(va[i], (double) rand());

					vb[i].mul_2si(vb[i], 53);
					vb[i].add(vb[i], (double) rand());
				}

			}


			// Chrono
			// ------

			fl_start = clock();

			g(tag, count, n[l], r, va, vb);

			fl_end = clock();


			fpnr_time += (difftime(fl_end, fl_start) / CLOCKS_PER_SEC);


			if ((K == 0) && (l == 0))
				cout << endl << "-------------------   " << tag << "   --------------------------" << endl;


			// Use va and r, see what has been involved in the computation
			for (int i = 0; i < n[l]; i++) {

				r.add(r, va[i]);

			}


			if (K == nbtrials - 1) {


				fpnr_time = fpnr_time / nbtrials;

				cout << r << endl; 


				double wflops = count * ((double) n[l]) / fpnr_time;
				cout << "flops     [" << n[l] << "]: " <<  wflops << endl;

				cout << "d ratio: " << dflops / wflops << endl << endl;

				if (l == 0)
					t1 = dflops / wflops;
				else if (l == nsize - 1)
					t2 = dflops / wflops;

			}

		} // Loop nbtrials

	} // Loop vector length

	cout << endl;
}


} // end namespace hplll

#endif



