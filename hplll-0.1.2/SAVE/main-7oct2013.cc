


#include "hlll.h"

using namespace hplll;

/****************************************************/
//
// Maquette de récursif 
//



  // U allouée à l'extérieur 

// ZT ici format pour les appels internes 
template<class ZT, class FT, class MatrixZT, class MatrixFT> 
void rec(ZZ_mat<mpz_t>& U, const ZZ_mat<mpz_t> B, 
	 Lattice<ZT, FT, MatrixZT, MatrixFT >& L, 	 
	 int l, int level, int lllmethod=HLLL) {

  int n,d,i,j;
  
  n=B.getRows();
  d=B.getCols();


  if (l <= 200) { // attention taille plus grande en fait pour double 

     if (lllmethod == HLLL) { 
       //cout << "Feuille " << endl; 
       L.put(B,1,0); 
       //print2maple(L.getbase(),n,d);
       L.hlll(0.99);
       U=L.getU();
     } 
     else if (lllmethod == FPLLL) { 
       
       ZZ_mat<mpz_t> AT;
       AT.resize(d,n);

       for (i=0; i<d; i++) 
	 for (j=0; j<n; j++) 
	   AT(i,j)=B(j,i);
       
      
       ZZ_mat<mpz_t> V;
       V.resize(d,d);
       for (i=0; i<d; i++) 
	 V(i,i)=1;

       //lllReduction(AT, 0.99, 0.51, LM_HEURISTIC,FT_DPE,0);
       lllReduction(AT, V, 0.99, 0.51, LM_WRAPPER,FT_DEFAULT,0);
      
       for (i=0; i<d; i++) 
	 for (j=0; j<d; j++) 
	   U(i,j)=V(j,i);

     } 
     //print2maple(U,d,d);

  } // Endif leaf 

   
   else {
     ZZ_mat<mpz_t> C1,C2;

     C1.resize(n,d);
     C2.resize(n,d);

     //cout << "Call A " << l/2 << endl; 
     trunc<mpz_t, ZZ_mat<mpz_t> >(C1, B, 1, d, n, l/2+d+2*(level+1), l/2); 
     //trunc<mpz_t, ZZ_mat<mpz_t> >(C1, B, 1, d, n, 0, l/2); 

     rec(U, C1, L, l/2, level+1, lllmethod);

     // ICI 
     if (level < 2) {
       cout << "*** " << maxbitsize(U) << "   " << maxbitsize(B) << endl; 
     }
     matprod(C1,B,U);
    

     ZZ_mat<mpz_t> U2;
     U2.resize(d,d);

     // cout << "Call B " << l/2 << endl; 
     trunc<mpz_t, ZZ_mat<mpz_t> >(C2, C1, 1, d, n, l/2+d+2*(level+1), 0); 
   
     rec(U2, C2, L, l/2, level+1, lllmethod);

     //rec(U2, C1, L, l/2, level+1, lllmethod);
     
     
     matprod(U,U2); 
   

   }
  

};


void recursive(ZZ_mat<mpz_t>& B, const ZZ_mat<mpz_t> A, int lllmethod=HLLL) { 


  //Lattice<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > > L(A,TRANSFORM,DEF_REDUCTION);
  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double,dpe_t> > L(A,TRANSFORM,DEF_REDUCTION);

   int n,d;
   int l;

   n=A.getRows();
   d=A.getCols();

   B.resize(n,d);

   ZZ_mat<mpz_t> U;
   U.resize(d,d);

   l=maxbitsize(A);

   rec(U, A, L, l, 0, lllmethod);

   matprod(B,A,U);

};


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

  /*  
  ZZ_mat<mpz_t> B;

  start = utime(); 
  startsec = utimesec(); 
  
  recursive(B, A, HLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  

  cout << endl; 

  if (output) print2maple(B,n,d);
 

  cout << endl; 
  cout << "--- cputime recursive hlll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime recursive hlll is " << tps_sec << " s" << endl;
  */
  /*
  start = utime(); 
  startsec = utimesec(); 
  
  recursive(B, A, FPLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  cout << endl; 
  cout << "--- cputime recursive fplll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime recursive fplll is " << tps_sec << " s" << endl;
  */

  double delta=0.9;

  // print2maple(A,n,d);
  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > L(A,NO_TRANSFORM,DEF_REDUCTION);

  start = utime(); 
  startsec = utimesec(); 
  
  L.hlll(delta);

  tps=utime()-start;
  tps_sec=utimesec()-startsec;
  
  cout << endl; 

  if (output) print2maple(L.getbase(),n,d);

  cout << endl; 
  cout << "--- cputime hlll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime hlll is " << tps_sec << " s" << endl;
  cout << "--- nblov = " << L.nblov <<  endl;
  cout << "--- compteur A = " << L.compteur <<  endl;
  cout << "--- compteur B = " << L.tmpcompt <<  endl;
  cout << "--- nbswaps  =  " << L.nbswaps <<  endl;

  //Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > >  L_check(L.getbase());
  //L_check.isreduced(delta);
  
  start = utime(); 
  startsec = utimesec(); 
  
  lllReduction(AT, delta, 0.51, LM_WRAPPER, FT_DEFAULT, 0);
  
  tps=utime()-start;
  tps_sec=utimesec()-startsec;
  
  cout << endl; 
  cout << "--- cputime fplll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime fplll is " << tps_sec << " s" << endl;

  if (output) {
    transpose(A,AT);
    print2maple(A,n,d);

  }

  return 0;
}
