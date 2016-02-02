


#include "../src/hlll.h"
//#include "matgen.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll;

#include "ratio.h"

int main(int argc, char *argv[])  {
  
  filebuf fb;
  iostream os(&fb);
  fb.open ("table.txt",ios::out);

  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------

  int k,N;

  int d=120;
  int bits=12000;
  fb.open ("table120.txt",ios::out);
   

  N=200;

  double delta=0.999;

  
  Timer time;

  os << endl << "LLL output quality test" << endl; 
  os <<         "-----------------------" << endl << endl;

  os << "With d = " << d << ",  bits = " << bits << ",  delta = " << delta <<  endl << endl;
  
   
  double  t,u,v,w,st=0.0,su=0.0,sv=0.0,sw=0.0,mt=0.0,mu=0.0,mv=0.0,mw=0.0;
  
  for (k=1; k<=N; k++) { 


      /*****************************************************************************/
      /*   k-th run  */
      /*****************************************************************************/
      
     

      A.resize(d+1,d);
      AT.resize(d,d+1);

      AT.gen_intrel(bits);
      transpose(A,AT);


      cout << "--------------  HLLL" << endl << endl; 

      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);  //* name 

      time.start();
      B.hlll(delta); //* name
      time.stop();

      ratio(B.getbase(),t,u,v,w);

      // Sums  
      st+=t;
      su+=u;
      sv+=v;
      sw+=w;

      mt=st/((double) k);
      mu=su/((double) k);
      mv=sv/((double) k);
      mw=sw/((double) k);

      os << endl << endl << "-------------- " << endl << endl;
      
      os << "Run " << k  << endl;

      os << endl << ".. log 2 Frobenius norm cond: " << t << endl;
      os << ".. Average diagonal ratio: " << u << endl;
      os << ".. Max diagonal ratio: " << v << endl;
      os << ".. First vector quality: " << w << endl;

      os << endl <<  "Current average " << endl;
      os << endl << ".. log 2 Frobenius norm cond: " << mt << endl;
      os << ".. Average diagonal ratio: " << mu << endl;
      os << ".. Max diagonal ratio: " << mv << endl;
      os << ".. First vector quality: " << mw << endl;
   
    }// End on runs, k loop


    // END 
    fb.close();


 
  return 0;
}
