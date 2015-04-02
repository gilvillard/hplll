#include "hplll/hlll.h"


using namespace std;
using namespace fplll;

using namespace hplll; 

int main(void) {

  typedef mpz_t integer_t;
  typedef matrix<Z_NR<integer_t> > MatrixZT;

  ZZ_mat<integer_t> A; // For hpLLL 
  ZZ_mat<integer_t> AT;  // fpLLL  

 
  A.resize(5,4); 
  AT.resize(4,5);  
  AT.gen_intrel(10);
  transpose(A,AT);

  Lattice<integer_t, double, MatrixZT, matrix<FP_NR<double> > > B(A,NO_TRANSFORM,DEF_REDUCTION);
  B.hlll(0.9);
  
  print2maple(B.getbase(),5,4);

}
