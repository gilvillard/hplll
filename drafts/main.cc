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


#include "hplll.h"

using namespace hplll;


int main(int argc, char *argv[])  {


	unsigned int old_cw;
	fpu_fix_start(&old_cw);

	int n = 8000;

	clock_t fl_start, fl_end;


	// DD
	// --

	double time;

	FP_NR<dd_real> f = 0.0;
	int prec = f.get_prec();

	vector<FP_NR<dd_real> > va;
	vector<FP_NR<dd_real> > vb;

	va.resize(n);
	vb.resize(n);


	for (int i = 0; i < n; i++) {

		va[i] = 0.0;
		vb[i] = 0.0;

		for (int k = prec; k > 0; k /= 53) {

			va[i].mul_2si(va[i], 53);
			va[i].add(va[i], (double) rand());

			vb[i].mul_2si(vb[i], 53);
			vb[i].add(vb[i], (double) rand());
		}

	}


	fl_start = clock();

	for (int i = 0; i < n; i++) {

		va[i].add(va[i], vb[i]);

	}

	fl_end = clock();


	// Use the result
	for (int i = 0; i < n; i++)
		f.add(f, va[i]);
	cout << f << endl;


	time = (difftime(fl_end, fl_start) / CLOCKS_PER_SEC);

	cout << endl << "Time dd: " << time << endl << endl;


	// QD
	// --

	double qtime;

	FP_NR<qd_real> qf = 0.0;
	int qprec = qf.get_prec();

	vector<FP_NR<qd_real> > qva;
	vector<FP_NR<qd_real> > qvb;

	qva.resize(n);
	qvb.resize(n);


	for (int i = 0; i < n; i++) {

		qva[i] = 0.0;
		qvb[i] = 0.0;

		for (int k = qprec; k > 0; k /= 53) {

			qva[i].mul_2si(qva[i], 53);
			qva[i].add(qva[i], (double) rand());

			qvb[i].mul_2si(qvb[i], 53);
			qvb[i].add(qvb[i], (double) rand());
		}

	}


	fl_start = clock();

	for (int i = 0; i < n; i++) {

		//qva[i].add(qva[i], qvb[i]);
		(qva[i]).get_data() = qd_real::sloppy_add((qva[i]).get_data(), (qvb[i]).get_data());

	}

	fl_end = clock();


	// Use the result
	for (int i = 0; i < n; i++)
		qf.add(qf, qva[i]);
	cout << qf << endl;


	qtime = (difftime(fl_end, fl_start) / CLOCKS_PER_SEC);

	cout << endl << "Time qd: " << qtime << endl << endl;

	cout << endl << "Ratio: " << qtime / time << endl << endl;


	return 0;
}








