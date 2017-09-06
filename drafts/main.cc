
#include <hplll.h>

using namespace hplll;

int main(int argc, char *argv[]) {


	Timer time, timed, timedd;

	ZZ_mat<mpz_t> A;

	ZZ_mat<mpz_t> AT;


	double delta = 0.99;

	int n, d;

	command_line_basis(A, n, d, delta, argc, argv);

	// fplll double

	AT.resize(d, n);
	transpose(AT, A);

	cout << AT << endl;

	timed.start();
	lll_reduction(AT, delta, 0.501, LM_FAST, FT_DOUBLE, 0, LLL_VERBOSE);
	timed.stop();


	// fplll dd

	AT.resize(d, n);
	transpose(AT, A);

	cout << AT << endl;

	timedd.start();
	lll_reduction(AT, delta, 0.501, LM_FAST, FT_DD, 0, LLL_VERBOSE);
	timedd.stop();


	// hplll

	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  B(A, NO_TRANSFORM, DEF_REDUCTION); //* name

	time.start();
	B.hlll(delta); //* name
	time.stop();


	cout << endl << endl << "   hplll double : " << time << endl << endl ;
	cout  << "   fplll double : " << timed << endl << endl ;
	cout << "   fplll dd : " << timedd << endl << endl ;
	cout << endl;











}
