
#include <hplll.h>


// *************************************************************
// requires output matrix C and vector delta to be allocated
// delta : ith denominators (minors), starting from 1
//
// todo: to d+1? 
// *************************************************************

int ff(ZZ_mat<mpz_t>& C, vector< Z_NR<mpz_t> >& delta, ZZ_mat<mpz_t> A) {

  int i,j,k;
  
  int d=A.getCols();

  matprod(C,transpose(A),A);

  print2maple(C,d,d);

  delta[0]=1;

  Z_NR<mpz_t> tz;
  
  for (i=0; i<d-1; i++) { // row that eliminates 

    for (k=i+1; k<d; k++) {

      // Row k with row i
      for (j=i+1; j<d; j++) {
	C(k,j).mul( C(i,i),C(k,j) );
	tz.mul( C(k,i), C(i,j) );
	C(k,j).sub(C(k,j),tz);

	mpz_divexact(C(k,j).GetData(),C(k,j).GetData(),delta[i].GetData());
	  
      }

      C(k,i)=0;
      
    }

    delta[i+1]=C(i,i);
    
  }

   print2maple(C,d,d);
   
  return 0;

}


// ***********************************************************
       
int main(int argc, char *argv[]) {
  
  ZZ_mat<mpz_t> A;

  int n,d;
  double delta_lll;

  command_line_basis(A, n, d, delta_lll, argc, argv); 

  cout << A << endl;

  ZZ_mat<mpz_t> C;
  C.resize(d,d);

  vector< Z_NR<mpz_t> > delta;
  delta.resize(d);

  ff(C,delta,A);

  // Unimodular 2x2 transformation 
  ZZ_mat<mpz_t> C;
  C.resize(d,d);


  
}
