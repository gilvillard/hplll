
#include <hplll.h>

#include "relations.h"

using namespace hplll;

int main(int argc, char *argv[]) {


	long alpha, d;

	filebuf fb;
	iostream os(&fb);

	vector<FP_NR<mpfr_t> > fpv;

	//  --------------   Test A

	// static string s;

	// fb.open ("C3_in", ios::in);

	// os >> alpha;
	// os >> d;

	// mpfr_set_default_prec(alpha);

	// fpv.resize(d);

	// for (int i = 0; i < d; i++) {
	// 	os >> s;
	// 	mpfr_set_str (fpv[i].get_data(), s.c_str(), 10, GMP_RNDN);
	// }

	// fb.close();

	// Z_mat<mpz_t> C;

	// FPTuple<long, double> L(fpv);

	// L.relation_f(C, alpha, 80, 20, 10, HLLL);

	// cout << C << endl;


	// ----------------  Test B

	
	// alpha = 2800;
	// mpfr_set_default_prec(alpha);

	// d = 65;
	// gen3r2s(fpv, d, 8, 8);



	// ZZ_mat<mpz_t> C;

	// FPTuple<long, double> L(fpv);

	// L.relation_f(C, alpha, 60, 200, 20, HLLL);
	// cout << C << endl;


	// ----------------  Test C

	
	alpha = 1800;
	mpfr_set_default_prec(alpha);

	d = 50;
	gen3r2s(fpv, d, 7, 7);



	ZZ_mat<mpz_t> C;

	FPTuple<mpz_t, dpe_t> L(fpv);

	L.relation_z(C, alpha,60,200,-1, FPLLL);
	cout << C << endl;











}






