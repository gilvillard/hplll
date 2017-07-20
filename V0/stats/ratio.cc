

#include "../src/hlll.h"

#include "ratio.h"


// ******************************************************************************
//
//  Main 
//
// ******************************************************************************

int main(void) {

  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT;  // fpLLL  

  /* For computing the gcd of 1136 and 2672 */

  int d=80;
  int bits = 1000;
  
  A.resize(d+1,d); 
  AT.resize(d,d+1);

  AT.gen_intrel(bits);
   
  transpose(A,AT);

  
  /* Via hlll lattice reduction */

  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >   L(A); 

  L.hlll(0.999);

  ZZ_mat<mpz_t> B;
  B.resize(d+1,d);
  B = L.getbase();
  
  // Computation of the ratio 

  double lfcond,av_ratio,max_ratio,c;
  
  print2maple(A,d+1,d);

  
  ratio(B,lfcond,av_ratio,max_ratio,c);


  cout << endl << endl << ".. log 2 Frobenius norm cond: " << lfcond << endl;
  cout << ".. Average diagonal ratio: " << av_ratio << endl;
  cout << ".. Max diagonal ratio: " << max_ratio << endl;
  cout << ".. First vector quality: " << c << endl; 
  //
  
  

}
