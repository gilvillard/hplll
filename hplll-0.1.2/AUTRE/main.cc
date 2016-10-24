

#include "fplll.h"


/* ***********************************************

          MAIN   

   ********************************************** */

using namespace fplll; 

int main(int argc, char *argv[])  {
  
  ZZ_mat<mpz_t>  A;

  A.read(cin);

  int start_time,run_time;

  start_time=cputime();
  
  lll_reduction(A, 0.99, 0.501, LM_FAST, FT_DEFAULT,0, LLL_VERBOSE);

  run_time=cputime()-start_time;


   cout << endl << "Time : " << run_time << endl << endl;

   
  return 0;
}
