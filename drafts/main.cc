
#include <hplll.h>

using namespace hplll;

int main(int argc, char *argv[]) {


	typedef __int128_t integer;

	Timer time, timed, timedd;

	ZZ_mat<mpz_t> A;


	double delta = 0.99;

	int n, d;

	//RandGen::init_with_time();

	command_line_basis(A, n, d, delta, argc, argv);


	ZZ_mat<integer> Along;

	matrix_cast(Along, A);

	// hplll

	//Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  B(A, NO_TRANSFORM, SEYSEN_REDUCTION); //* name
	Lattice < integer, double, matrix<Z_NR<integer> >, matrix<FP_NR<double> > >  B(Along, NO_TRANSFORM, SEYSEN_REDUCTION); //* name

	verboseDepth = 0;

	time.start();
	B.hlll(delta); //* name
	time.stop();

	// RE-USE OF A
	matrix_cast(A, B.getbase());

	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A); //* names

	T.isreduced(delta - 0.1);

	double t, u, v, w;

	hplll::ratio<mpz_t>(A, t, u, v, w);

	cout << endl << ".. fplll log 2 Frobenius norm cond: " << t << endl;
	cout << ".. Average diagonal ratio: " << u << endl;
	cout << ".. Max diagonal ratio: " << v << endl;
	cout << ".. First vector quality: " << w << endl;


	cout << endl << endl << "   hplll double : " << time << endl ;



	cout << endl;











}
