
#include <hplll.h>


using namespace hplll;

int main(int argc, char *argv[]) {


	long alpha, d;

	filebuf fb;
	iostream os(&fb);

	vector<FP_NR<mpfr_t> > fpv;

	static string s;

	cin  >> alpha;
	cin  >> d;

	mpfr_set_default_prec(alpha);

	fpv.resize(d);

	for (int i = 0; i < d; i++) {
		cin >> s;
		mpfr_set_str (fpv[i].get_data(), s.c_str(), 10, GMP_RNDN);
	}

	ZZ_mat<mpz_t> C;

	FPTuple<long, double, matrix<FP_NR<double> > > L(fpv);

	L.relation(C, alpha, 10, 10, 40, FPLLL);

	cout << C << endl;

}






