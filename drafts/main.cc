
#include <hplll.h>

using namespace hplll;




void go(int n, int K, double alpha) {


	RandGen::init_with_time();

	vector<double> lcond(K);

	double scond = 0.0;

	double t, u, v, w;

	double delta = 0.99;

	Timer time;

	int status;

	typedef __int128_t integer;
	//typedef long integer;

	// Storing anomalies

	int ano = 0;

	char* s;

	s = (char *) malloc(600);

	filebuf fb;
	iostream os(&fb);


	cout << endl << endl << "*************************    n = " << n << "   alpha = " << alpha << endl;

	ZZ_mat<mpz_t> A;
	A.resize(n, n);

	ZZ_mat<integer> Along;



	cout << endl;

	for (int k = 0; k < K; k++) {

		genalpha(A, n, alpha);

		matrix_cast(Along, A);



		Lattice < integer, double, matrix<Z_NR<integer> >, matrix<FP_NR<double> > >  B(Along, NO_TRANSFORM, SEYSEN_REDUCTION);

		verboseDepth = 1;

		time.clear();
		time.start();

		B.hlll(delta);

		time.stop();

		verboseDepth = 0;

		// RE-USE OF A
		matrix_cast(A, B.getbase());

		// Check only for the first trial
		// ------------------------------

		hplll::ratio<mpz_t>(A, t, u, v, w);

		//if ((k == 0) || (k == K - 1))
		{


			Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A);

			status = T.isreduced(delta - 0.1);


			if (status == -1) {

				strncpy(s, "ano.txt", 599);

				sprintf(s + 4, "%d", ano);
				ano += 1;

				fb.open (s, ios::out);
				os << "Dimension: " << n << endl;
				os << A;
				fb.close();

				print2maple(A, n, n);

			}

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


	cout << endl << "Average log 2 Frobenius norm cond : " << scond / K << endl;

	double ln;
	ln = log(n) / log(2.0);
	ln *= ln;

	cout << endl << "Average over log 2 n^2 : " << (scond / K) / ln << endl;


	cout << endl << "Time: " << time << endl ;

}


int main(int argc, char *argv[]) {

	int n = 420;

	double delta = 0.99;

	
	typedef __int128_t integer;
	//typedef long integer;

	ZZ_mat<mpz_t> A;
	A.resize(n, n);

	cin >> A;


	ZZ_mat<integer> Along;

	matrix_cast(Along, A);


	Lattice < integer, double, matrix<Z_NR<integer> >, matrix<FP_NR<double> > >  B(Along, NO_TRANSFORM, SEYSEN_REDUCTION);

	verboseDepth = 1;

	B.hlll(delta);

	matrix_cast(A, B.getbase());



}








