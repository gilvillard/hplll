


#include <hplll.h>


using namespace hplll;

int main(int argc, char *argv[]) {

	Timer time;

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


	// alpha = 1600;
	// mpfr_set_default_prec(alpha);

	// int r, s;
	// r = 7;
	// s = r;

	// d = r*s+1;
	// gen3r2s(fpv, d, r, s);



	// ZZ_mat<mpz_t> C;



	// // cout << alpha << endl;
	// // cout << d << endl;
	// // for (int i = 0; i < d; i++) {
	// // 	//mpfr_out_str (stdout, 10, alpha, fpv[i].get_data(), GMP_RNDN);
	// // 	mpfr_printf ("%.1940Rf", fpv[i].get_data());
	// // 	cout << endl;
	// // }


	// FPTuple<mpz_t, dpe_t, MatrixPE<double, dpe_t> > L(fpv);

	// time.start();

	// L.lll(C, alpha);

	// //L.relation_f(C, alpha, 60, 100, 20, FPLLL);

	// time.stop();

	// cout << C << endl;

	// cout << endl << endl << "   relation : " << time << endl ;

	//  --------------   Test from file / Poisson

	static string s;



	//fb.open ("alpha.in", ios::in);

	//os >> alpha;
	//os >> d;

	cin  >> alpha;
	cin  >> d;

	mpfr_set_default_prec(alpha);

	fpv.resize(d);

	for (int i = 0; i < d; i++) {
		//os >> s;
		cin >> s;
		mpfr_set_str (fpv[i].get_data(), s.c_str(), 10, GMP_RNDN);
	}

	//fb.close();

	ZZ_mat<mpz_t> C;


	FPTuple<__int128_t, double, matrix<FP_NR<double> > > L(fpv);




#ifdef _OPENMP
	double st = omp_get_wtime();
	OMPTimer ptime;
	ptime.start();
#endif 

	time.start();


	L.set_num_threads(2);

	L.relation(C, alpha, 20, 60, 100, HLLL);


#ifdef _OPENMP
	double en = omp_get_wtime();
	ptime.stop();
#endif 
	time.stop();



	cout << C << endl;

#ifdef _OPENMP
	cout << endl << endl << "   omp time relation : " << en - st  << endl ;
	cout << endl << "   omp timer relation : " << ptime << endl ;
#endif 

	cout << endl << "   timer relation : " << time << endl << endl;


}






