

#include "fplll.h"


/* ***********************************************

          MAIN   

   ********************************************** */

int main(int argc, char *argv[])  {


  int start_time,run_time;

   
   ZZ_mat<mpz_t>  L;

   L.read(cin);

   start_time=cputime();
  
   lllReduction(L, 0.99, 0.501, LM_FAST, FT_DEFAULT,0, LLL_VERBOSE);

   run_time=cputime()-start_time;


   cout << endl << "Time : " << run_time << endl << endl;

  
   
  return 0;
}
