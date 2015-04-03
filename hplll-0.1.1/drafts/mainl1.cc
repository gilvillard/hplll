


#include "hlll.h"
#include "l1.h" 

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

    n=atoi(argv[decal+1]);
    d=n-1;

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

  
  ZZ_mat<mpz_t> B;

  start = utime(); 
  startsec = utimesec(); 
  
  l1(B, A, 800, HLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  

  cout << endl; 

  if (output) print2maple(B,n,d);
 

  cout << endl; 
  cout << "--- cputime l1 hlll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime l1 hlll is " << tps_sec << " s" << endl;
  
  
  
  start = utime(); 
  startsec = utimesec(); 
  
  l1(B, A, 800, FPLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  cout << endl; 
  cout << "--- cputime l1 fplll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime l1 fplll is " << tps_sec << " s" << endl;
  
  
  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L(A,NO_TRANSFORM,DEF_REDUCTION);

  start = utime(); 
  startsec = utimesec(); 
  
  L.hlll(0.99);

  tps=utime()-start;
  tps_sec=utimesec()-startsec;
  
  cout << endl; 

  if (output) print2maple(L.getbase(),n,d);

  cout << endl; 
  cout << "--- cputime hlll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime hlll is " << tps_sec << " s" << endl;
   
  //Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > >  L_check(L.getbase());
  //L_check.isreduced(0.8);
  
  
  start = utime(); 
  startsec = utimesec(); 
  
  lllReduction(AT, 0.99, 0.51, LM_WRAPPER, FT_DEFAULT, 0);
  
  tps=utime()-start;
  tps_sec=utimesec()-startsec;
  
  cout << endl; 
  cout << "--- cputime fplll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime fplll is " << tps_sec << " s" << endl;
  


  return 0;
}
