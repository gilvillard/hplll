

#include "hlll.h"
#include "ratio.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  double t,u,v,w;
  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------
  { 
  
    cout << "************************************************************************** " << endl; 
    int d=80;
    int nbbits=8000;

    Timer time;

    double delta=0.999;

   
    A.resize(d+1,d); 
    tmpmat.resize(d+1,d); 
    AT.resize(d,d+1);  
    AT.gen_intrel(nbbits);
    transpose(A,AT);

    cout << "--------------  SeysenLLL" << endl << endl; 

   
    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,SEYSEN_REDUCTION);

    time.start();
    B.hlll(delta);
    time.stop();
    
    ratio(B.getbase(),t,u,v,w);
    
    cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
    cout << ".. Average diagonal ratio: " << u << endl;
    cout << ".. Max diagonal ratio: " << v << endl;
    cout << ".. First vector quality: " << w << endl;
   
    cout << endl << "Time: " << time << endl;    

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T1.isreduced(delta);

    cout << endl; 

    
    cout << "--------------  HLLL" << endl << endl; 
   

    Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > C(A,NO_TRANSFORM,DEF_REDUCTION);
    time.start();
    C.hlll(delta);
    time.stop();

    //print2maple(C.getbase(),d+1,d);

    ratio(C.getbase(),t,u,v,w);
    
    cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
    cout << ".. Average diagonal ratio: " << u << endl;
    cout << ".. Max diagonal ratio: " << v << endl;
    cout << ".. First vector quality: " << w << endl;

    cout << endl << "Time: " << time << endl;
     
    //cout << "   bits = " << nbbits << endl;
    //cout << "   dimension = " << d  << endl;
    //cout << "   time B: " << start/1000 << " ms" << endl;
    //cout << "   time B: " << startsec << " s" << endl;


    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(C.getbase(),NO_TRANSFORM,DEF_REDUCTION);
    T2.isreduced(delta);

    /*
    
    cout << endl; 
  
    transpose(AT,A);
 
    start=utime();
    startsec=utimesec();
    lllReduction(AT, delta, 0.5, LM_FAST,FT_DEFAULT,0);
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    cout << "   bits = " << nbbits << endl;
    cout << "   dimension = " << d  << endl;
    cout << "   time C: " << start/1000 << " ms" << endl;
    cout << "   time C: " << startsec << " s" << endl;
  
    transpose(tmpmat,AT);
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T3(tmpmat,NO_TRANSFORM,DEF_REDUCTION);
    T3.isreduced(delta);
    */

  } 

 
  return 0;
}
