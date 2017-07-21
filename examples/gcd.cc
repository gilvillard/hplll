
#include <hplll.h>


int main(void) {

  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT;  // fpLLL  

  /* For computing the gcd of 1136 and 2672 */
  
  A.resize(2,2); 
  AT.resize(2,2);

  A(0,0)=113600;
  A(0,1)=267200;
  A(1,0)=1;
  A(1,1)=0;
  
  transpose(AT,A);

  /* Via hlll lattice reduction */
  
  Lattice<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > > B(A);
  B.hlll(0.99);

  
  print2maple(B.getbase(),2,2);

  cout << endl << B.getbase() << endl; 

  lll_reduction(AT, 0.99, 0.501, LM_WRAPPER);

  cout << endl << AT << endl; 
  

}
