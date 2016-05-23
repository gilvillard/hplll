

#include "hplll.h"


/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll;

int main(int argc, char *argv[]) {

  //typedef mpz_t ZT;
  //typedef double ZT;

  int d=20;
  int nbbits=60;

  int K=200;
  
  // PARSE_MAIN_ARGS {
    
  //   MATCH_MAIN_ARGID("-d",d);
  //   MATCH_MAIN_ARGID("-bits",nbbits);
  // }

  for (d=2; d<10000; d++) {

    OMPTimer tots,totp;
    tots.clear();
    totp.clear();

    for (int k=0; k<K; k++) {
    
      Matrix<FP_NR<double> > A,B,C,D;
  
      A.resize(d,d); 
      B.resize(d,d);
      C.resize(d,d);
      D.resize(d,d);

      Z_NR<mpz_t> tz;
      
      for (int i=0; i<d; i++)
	for (int j=0; j<d ; j++) {
	  tz.randb(nbbits);
	  A(i,j).set_z(tz);
	}

      

      for (int i=0; i<d; i++)
	for (int j=0; j<d ; j++) {
	  tz.randb(nbbits);
	  B(i,j).set_z(tz);
	}
  

      OMPTimer tseq,tpar;

      tseq.clear();
      tpar.clear();

      // SEQ
      // ---


      tseq.start();

    
      matprod(C,A,A);
      matprod(D,B,B);
    
 
      tseq.stop();
      tots+=tseq;
      
      // PAR
      // ---

      tpar.start();
  
#pragma omp parallel sections num_threads(2)
      {
	
#pragma omp section
	{
	  matprod(C,A,A);
	}
    
#pragma omp section
	{
	  matprod(D,B,B);
	}
      }
  
      tpar.stop();
      totp+=tpar;
      

    }
    
    cout << endl << "-------- Seq --------" << endl << endl;

    

    cout << endl << tots << endl;

    cout << endl << "-------- Par --------" << endl << endl;

    
    cout << endl << totp << endl;


    double ratio= totp.realtime()/tots.realtime();
    
    cout << endl << endl << "Dim: " << d << "     Ratio par/seq: " << ratio << endl << endl; 
    cout << "               Par: " <<  totp.realtime()/((double) K) << endl;
    cout << "               Seq: " <<  tots.realtime()/((double) K) << endl;
  
    if (ratio < 0.5601) d=2000000; 
  }
  
  
}


