


#include "hlll.h"
#include "l1.h" 
#include "block.h" 

using namespace hplll;


/* --------------------------------------------------*/

int main(int argc, char *argv[])  {
 

  int n=0,d=0; 
  int decal=1,setprec;

  /* ----------------------------------------- */
  typedef mpz_t integer_t;

  
  ZZ_mat<integer_t> A; 
  ZZ_mat<integer_t> AT;  

  int output=0;

 
  if (strcmp(argv[decal],"-k")==0) {

    d=atoi(argv[decal+1]);
    n=d+1;

    A.resize(n,d); 
    AT.resize(d,n);  

    AT.gen_intrel(atoi(argv[decal+2]));
    transpose(A,AT);
    output = atoi(argv[decal+3]);

  }

  if (argc >= 6) { 
    if (strcmp(argv[decal+4],"-prec")==0) {
  
      setprec = atoi(argv[decal+5]);
      
    }
  };

  // ---------------------------------------------------
  
  int start,startsec,tps=0,tps_sec=0;
  

  if (output) print2maple(A,n,d);

  //------------------
  ZZ_mat<mpz_t> B1;

  start = utime(); 
  startsec = utimesec(); 

  B1.resize(n,d);  
  blevel(B1, A, 2, HLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  cout << "--- cputime blll 1 is " << tps/1000 << " ms" << endl;
  cout << "--- cputime blll 1 is " << tps_sec << " s" << endl;

  cout << endl; 

  if (output) print2maple(B1,n,d);
  //---------------------
  ZZ_mat<mpz_t> B2;

  start = utime(); 
  startsec = utimesec(); 

  B2.resize(n,d);  
  blevel(B2, B1, 2, HLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  cout << "--- cputime blll 2 is " << tps/1000 << " ms" << endl;
  cout << "--- cputime blll 2 is " << tps_sec << " s" << endl;


  cout << endl; 

  if (output) print2maple(B2,n,d);
  //---------------------
  ZZ_mat<mpz_t> B3;

  start = utime(); 
  startsec = utimesec(); 

  B3.resize(n,d);  
  blevel(B3, B2, d, HLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  cout << "--- cputime blll 3 is " << tps/1000 << " ms" << endl;
  cout << "--- cputime blll 3 is " << tps_sec << " s" << endl;

  cout << endl; 

  if (output) print2maple(B3,n,d);
  //---------------------

 //---------------------

  start = utime(); 
  startsec = utimesec(); 


  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(B2,NO_TRANSFORM,DEF_REDUCTION);
  C.hlll(0.99);

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  cout << "--- cputime blll 3b is " << tps/1000 << " ms" << endl;
  cout << "--- cputime blll 3b is " << tps_sec << " s" << endl;

  cout << endl; 

 
  cout << endl; 
  
  
  /* L1 FPLLL 

  start = utime(); 
  startsec = utimesec(); 
  
  l1(B, A, 800, FPLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  cout << endl; 
  cout << "--- cputime l1 fplll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime l1 fplll is " << tps_sec << " s" << endl;
  */

  /* HLLL */ 
  
  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L(A,NO_TRANSFORM,DEF_REDUCTION);

  start = utime(); 
  startsec = utimesec(); 
  
  L.hlll(0.99);

  tps=utime()-start;
  tps_sec=utimesec()-startsec;
  
  cout << endl; 

  if (output) print2maple(L.getbase(),n,d);

  cout << endl; 
  cout << "nblov: " << L.nblov << endl; 
  cout << "--- cputime hlll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime hlll is " << tps_sec << " s" << endl;
   
  //Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > >  L_check(L.getbase());
  //L_check.isreduced(0.8);
  

  /* FPLLL 
  
  start = utime(); 
  startsec = utimesec(); 
  
  lllReduction(AT, 0.99, 0.51, LM_WRAPPER, FT_DEFAULT, 0);
  
  tps=utime()-start;
  tps_sec=utimesec()-startsec;
  
  cout << endl; 
  cout << "--- cputime fplll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime fplll is " << tps_sec << " s" << endl;
  
  */


  return 0;
}
