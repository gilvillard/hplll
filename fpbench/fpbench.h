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

#include <chrono>

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

	for (int l = 0; l < nsize; l++) {

		va.resize(n[l]);
		vb.resize(n[l]);


		Z_NR<mpz_t> r0;

		srand(time(NULL));
		r0.randb(prec);




		string tag;

		Timer t;
		clock_t fl_start, fl_end;

		//int nbtrials = 20;

		double fpnr_time[nbtrials], fp_time[nbtrials];

		auto fpnrbegin = chrono::high_resolution_clock::now();
		auto fpnrend =  chrono::high_resolution_clock::now();
		auto fpnrduration = 0.0;

		auto fpbegin = chrono::high_resolution_clock::now();
		auto fpend = chrono::high_resolution_clock::now();
		auto fpduration = 0.0;



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


			t.clear();

			t.start();
			fl_start = clock();

			fpnrbegin = chrono::high_resolution_clock::now();

			g(tag, count, n[l], r, va, vb);

			fpnrend = chrono::high_resolution_clock::now();

			fl_end = clock();
			t.stop();

			fpnrduration += chrono::duration_cast<chrono::nanoseconds>(fpnrend - fpnrbegin).count();



			if ((K == 0) && (l == 0))
				cout << endl << "-------------------   " << tag << "   --------------------------" << endl;

			fpnr_time[K] = difftime(fl_end, fl_start) / CLOCKS_PER_SEC;

			// Use va and r, see what has been involved in the computation
			for (int i = 0; i < n[l]; i++) {

				r.add(r, va[i]);

			}

			// Pure double speed

			vector<double> r1;
			vector<double> r2;

			r1.resize(n[l]);
			r2.resize(n[l]);

			for (int i = 0; i < n[l]; i++) {

				r1[i] = rand();
				r2[i] = rand();

			}

			fl_start = clock();

			fpbegin = chrono::high_resolution_clock::now();

			for (int i = 0; i < n[l]; i++) {

				r1[i] += 3416 * r1[i]; // * r2[i];

			}

			fpend = chrono::high_resolution_clock::now();

			fl_end = clock();


			fpduration += chrono::duration_cast<chrono::nanoseconds>(fpend - fpbegin).count();


			for (int i = 0; i < n[l]; i++)
				fpr += r1[i];

			fp_time[K] = difftime(fl_end, fl_start) / CLOCKS_PER_SEC;



			if (K == nbtrials - 1) {

				double avg_fpnr = 0.0;
				double avg_fp = 0.0;

				for (int K = 0; K < nbtrials; K++) {

					avg_fpnr += fpnr_time[K];
					avg_fp += fp_time[K];

				}

				avg_fpnr /= nbtrials;
				avg_fp /= nbtrials;

				cout << fpr << endl;
				cout << r << endl;


				//cout << "CPU & clock time:  " << t.usertime() << "  "  << avg_fpnr  << endl << endl;
				//cout << "CPU & clock time:  " << t.usertime() << "  "  << avg_fp << endl << endl;

				//cout << "flops     [" << n[l] << "]: " << count * ((double) n[l]) / t.usertime() << endl;

				cout << endl << tag << endl << endl;

				cout << fpnrduration << "ns total, average : " << fpnrduration / nbtrials << "ns." << endl;

				fpnrduration /= nbtrials;
				fpduration /= nbtrials;

				cout << "GFlops     [" << n[l] << "]: " << count * ((double) n[l]) / fpnrduration << endl;
				cout << "flops     [" << n[l] << "]: " << count * ((double) n[l]) / avg_fpnr << endl;
				cout << "dGFlops    [" << n[l] << "]: " << 2 * ((double) n[l]) / fpduration << endl;
				//cout << "dflops    [" << n[l] << "]: " << 2 * ((double) n[l]) / avg_fp << endl;

				//cout << "d ratio: " << 2 * avg_fpnr / (count * avg_fp) << endl << endl;
				cout << "d ratio: " << ((double) 2) * fpnrduration / (((double) count) * fpduration) << endl << endl;

				if (l == 0)
					t1 = ((double) 2) * fpnrduration / (((double) count) * fpduration);
				else if (l == nsize - 1)
					t2 = ((double) 2) * fpnrduration / (((double) count) * fpduration);

			}

		} // Loop nbtrials

	} // Loop vector length

	cout << endl;
}


} // end namespace hplll

#endif


