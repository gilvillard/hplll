#include "../src/hlll.h"


using namespace std;
using namespace fplll;

using namespace hplll; 

// ******************************************************************************
//
//  Statitics concerning the quality of the reduction
//  The matrix is assumed to be reduced   n rows >= d columns
//
//   log_2 of the Frobenius norm cond 
//
// ******************************************************************************


int ratio(ZZ_mat<mpz_t> B, double& lfcond) {

  int i,j,k;
  
  cout << endl << B << endl;

  //int n=B.getRows();
  int d=B.getCols();
  
  mpfr_set_default_prec(2*d);
  
  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(B);

  L.householder();

  matrix<FP_NR<mpfr_t> > R;

  R= L.getR();

  // ** Frobenius norm cond  ** 
  // **************************

  matrix<FP_NR<mpfr_t> > aR;
  matrix<FP_NR<mpfr_t> > iR;

  // Absolute value
  
  aR.resize(d,d);
      
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      aR(i,j).abs(R(i,j));

  // Inversion 
  
  iR.resize(d,d);
  
  for (i=0; i<d; i++) iR(i,i)=1.0; 
  
  // Higham p263 Method 1
  
  FP_NR<mpfr_t> t,one;
  one=1.0;

  for (j=0; j<d; j++) {
    iR(j,j).div(one,aR(j,j));
    for (k=j+1; k<d; k++) {
      t.mul(iR(j,j),R(j,k));
      iR(j,k).neg(t);
    }
    
    for (k=j+1; k<d; k++) {
      for (i=j+1; i<k; i++) {
	t.mul(iR(j,i),R(i,k)); 
	iR(j,k).sub(iR(j,k),t); 
      }
      iR(j,k).div(iR(j,k),R(k,k)); 
    }
  }

  // Abs iR --> iR 
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      iR(i,j).abs(iR(i,j));
  
  // aR * iR 
  
  matrix<FP_NR<mpfr_t> > prod;
  prod.resize(d,d);
  
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) {
      prod(i,j).mul(aR(i,0),iR(0,j));
      for (k=1; k<d ; k++) 
	prod(i,j).addmul(aR(i,k),iR(k,j));
    }


  // Frobenius norm 

  FP_NR<mpfr_t>  cc;
  cc=0.0; 
  
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      cc.addmul(prod(i,j),prod(i,j));
  
  cc.sqrt(cc); 

  mpfr_log2(cc.getData(),cc.getData(),GMP_RNDN);

  lfcond = cc.get_d();
  
  

  print2maple(R,d,d);
    
  return 0;
  
}

// ******************************************************************************
//
//  Main 
//
// ******************************************************************************

int main(void) {

  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT;  // fpLLL  

  /* For computing the gcd of 1136 and 2672 */

  int d=8;
  int bits = 20;
  
  A.resize(d+1,d); 
  AT.resize(d,d+1);

  AT.gen_intrel(bits);
   
  transpose(A,AT);

  
  /* Via hlll lattice reduction */
  
  Lattice<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > > L(A);
  L.hlll(0.99);

  ZZ_mat<mpz_t> B;
  B.resize(d+1,d);
  B = L.getbase();
  
  // Computation of the ratio 

  double lfcond;
  
  print2maple(A,d+1,d);

  
  ratio(B,lfcond);


  cout << endl << endl << ".. log 2 Frobenius norm cond: " << lfcond << endl; 
  //
  
  

}
