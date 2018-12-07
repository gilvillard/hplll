
#include <hplll.h>
#include <slll.h>

using namespace hplll;

int main(int argc, char *argv[]) {


	Timer time, timep, timed, timedd;

	ZZ_mat<mpz_t> A;

	ZZ_mat<mpz_t> AT;


	double delta = 0.99;

	int n, d;

	command_line_basis(A, n, d, delta, argc, argv);

	// fplll double

	AT.resize(d, n);
	transpose(AT, A);

	//cout << AT << endl;

	timed.start();
	//lll_reduction(AT, delta, 0.501, LM_FAST, FT_DOUBLE, 0, LLL_VERBOSE);
	timed.stop();


	// // fplll dd

	// AT.resize(d, n);
	// transpose(AT, A);

	// //cout << AT << endl;

	// timedd.start();
	// lll_reduction(AT, delta, 0.501, LM_FAST, FT_DD, 0, LLL_VERBOSE);
	// timedd.stop();


	// hplll

	//SLattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  B(A, 4); //* name

	int S = 32;

	SLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A, S);

	B.setprec(1600);

	timep.start();

	B.phouseholder(S);

	timep.stop();

	time.start();

	B.phouseholder(1);

	time.stop();

	cout << endl << endl << "   householder p : " << timep << endl;
	cout << endl << "   householder : " << time << endl;
	cout << endl << "   ratio : " << ((double) time.realtime()) / ((double) timep.realtime()) << endl << endl;

	//verboseDepth = 1;

	// time.start();
	// B.hlll(delta); //* name
	// time.stop();

	// transpose(A, AT);
	// Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(B.getbase()); //* names

	// T.isreduced(delta - 0.1);

	// double t, u, v, w;

	// hplll::ratio<mpz_t>(B.getbase(), t, u, v, w);

	// cout << endl << ".. fplll log 2 Frobenius norm cond: " << t << endl;
	// cout << ".. Average diagonal ratio: " << u << endl;
	// cout << ".. Max diagonal ratio: " << v << endl;
	// cout << ".. First vector quality: " << w << endl;


	// cout << endl << endl << "   hplll double : " << time << endl ;
	// cout << "   Householder: " << B.dbg << endl << endl ;
	// cout  << "   fplll double : " << timed << endl << endl ;
	// //cout << "   fplll dd : " << timedd << endl << endl ;
	// cout << endl;











}
