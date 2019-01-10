
#include <hplll.h>

using namespace hplll;




void go(int n, int K, double alpha) {


	RandGen::init_with_time();

	vector<double> lcond(K);

	double scond = 0.0;

	double t, u, v, w;

	double delta = 0.99;

	Timer time;


	typedef __int128_t integer;
	//typedef long integer;


	cout << endl << endl << "*************************    n = " << n << "   alpha = " << alpha << endl;

	ZZ_mat<mpz_t> A;
	A.resize(n, n);

	ZZ_mat<integer> Along;


	time.start();


	cout << endl;

	for (int k = 0; k < K; k++) {

		genalpha(A, n, alpha);

		matrix_cast(Along, A);

		Lattice < integer, double, matrix<Z_NR<integer> >, matrix<FP_NR<double> > >  B(Along, NO_TRANSFORM, SEYSEN_REDUCTION);

		verboseDepth = 0;

		B.hlll(delta);

		// RE-USE OF A
		matrix_cast(A, B.getbase());

		// Check only for the first trial
		// ------------------------------

		hplll::ratio<mpz_t>(A, t, u, v, w);

		if ((k == 0) || (k==K-1)) {

		
			Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A);
			T.isreduced(delta - 0.1);

			cout << endl << ".. fplll log 2 Frobenius norm cond: " << t << endl;
			cout << ".. Average diagonal ratio: " << u << endl;
			cout << ".. Max diagonal ratio: " << v << endl;
			cout << ".. First vector quality: " << w << endl;
			cout << ".. Nb swaps: " << B.nbswaps << endl;

		}


		lcond[k] = t;

		scond += t;

		cout.flush() << t << ",";

	}
	cout << endl;

	time.stop();

	cout << endl << "Average log 2 Frobenius norm cond : " << scond / K << endl;

	double ln;
	ln = log(n) / log(2.0);
	ln *= ln;

	cout << endl << "Average over log 2 n^2 : " << (scond / K) / ln << endl;


	cout << endl << "Time: " << time << endl ;

}


int main(int argc, char *argv[]) {

	int n; // Initial value

	for (int i = 0; i < 100; i++ ) {

		n = 20 + i * 40;

		go(n, 20, 1.04);

	}


}








