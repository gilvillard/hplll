
#include <hplll.h>

using namespace hplll;

int main(int argc, char *argv[]) {


	Timer time, timed, timedd;

	ZZ_mat<mpz_t> A0;
	ZZ_mat<mpz_t> A;

	ZZ_mat<mpz_t> AT;


	double delta = 0.99;

	int n, d;

	command_line_basis(A0, n, d, delta, argc, argv);

	// fplll double

	A.resize(d, d);
	for (int i = 0; i < d; i++)
		for (int j = 0; j < d; j++)
			A(i, j) = A0(i, j);

	AT.resize(d, d);
	transpose(AT, A);

	//cout << AT << endl;


	// Write the inpute matrix in a file
	filebuf fb;
	iostream os(&fb);
	fb.open ("/Users/gvillard/Installed-Libraries/code/src/lat.txt", ios::out);
	os << AT;
	fb.close();


	timed.start();
	lll_reduction(AT, delta, 0.501, LM_FAST, FT_DOUBLE, 0, LLL_VERBOSE);
	//lll_reduction(AT, delta, 0.501, LM_FAST, FT_DOUBLE, 0);
	timed.stop();


	// // fplll dd

	// AT.resize(d, n);
	// transpose(AT, A);

	// //cout << AT << endl;

	// timedd.start();
	// lll_reduction(AT, delta, 0.501, LM_FAST, FT_DD, 0, LLL_VERBOSE);
	// timedd.stop();


	// hplll

	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  B(A, NO_TRANSFORM, DEF_REDUCTION); //* name

	verboseDepth = 1;

	time.start();
	B.hlll(delta); //* name
	time.stop();

	cout << endl << endl;

	verboseDepth = 0;

	transpose(A, AT);
	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(B.getbase()); //* names

	T.isreduced(delta - 0.1);

	double t, u, v, w;

	hplll::ratio<mpz_t>(B.getbase(), t, u, v, w);

	cout << endl << ".. fplll log 2 Frobenius norm cond: " << t << endl;
	cout << ".. Average diagonal ratio: " << u << endl;
	cout << ".. Max diagonal ratio: " << v << endl;
	cout << ".. First vector quality: " << w << endl;


	cout << endl << endl << "   hplll double : " << time << endl ;
	cout << "   Householder: " << B.dbg << endl << endl ;
	cout  << "   fplll double : " << timed << endl << endl ;
	//cout << "   fplll dd : " << timedd << endl << endl ;
	cout << endl;



}
