

#include "hlll.h"

#include <omp.h>



/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll;

int main() {

  int d=32;  // Divisible par S 
  int S=4;
  int Sdim=d/S;

  int nbbits=10;
  int i,j,k;


  matrix<Z_NR<mpz_t> > B;
  B.resize(d,d);

  matrix<Z_NR<mpz_t> > C;
  C.resize(d,d);
 

  Matrix<Z_NR<mpz_t> > resB,resC;
  resB.resize(d,d);
  resC.resize(d,d);

  Matrix<Z_NR<mpz_t> > U;
  U.resize(d,d);

  for (i=0; i<d; i++) 
    for (j=0; j<d; j++) 
      (B(i,j)).randb(nbbits);
  
  set(C,B);  


  for (k=0; k<S ; k++) {

    int dec=k*Sdim;

    for (i=0; i<Sdim; i++) 
      for (j=0; j<d-dec; j++) 
	(U(dec+i,dec+j)).randb(nbbits);

  }


  // les deux calculs 

  print2maple(B,d,d);
 
  print2maple(U,d,d);

  pmatprod_in(B,U,S);

  pmatprod(C,U,S,3*S); 


  set(resB,B);
  set(resC,C);
  

  // Comparaison 

  if (matcmp(resB,resC,d,d) ==1) cout << endl << "***** Ok " << endl; 
  else cout << endl << "***** Error" << endl; 
  
}


