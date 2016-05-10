

#include "hlll.h"

#include <omp.h>



/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll;

int main() {

  typedef long ZT;
  
  int d=400;
  int nbbits=32;
  int i;

  OMPTimer time,ptime;
  
  vector<Z_NR<ZT> > u,v,w,z;
  Z_NR<ZT> alpha,beta;
  
  u.resize(d);
  v.resize(d);
   
  for (i=0; i<d; i++) u[i].randb(nbbits);
  for (i=0; i<d; i++) v[i].randb(nbbits); 
  
  w.resize(d);
  z.resize(d);
   
  for (i=0; i<d; i++) w[i].randb(nbbits);
  for (i=0; i<d; i++) z[i].randb(nbbits); 

  
  alpha.randb(20);
  beta.randb(20);

   // ************************
  
   time.start();

   for (i=0; i<d; i++)
     v[i].addmul(u[i],alpha);

   for (i=0; i<d; i++)
     z[i].addmul(w[i],beta);

  
   time.stop();
   
   cout << endl << "-------- Seq --------" << endl << endl;

   cout << endl << time << endl;

   // ************************

   u.resize(d);
   v.resize(d);
   
   for (i=0; i<d; i++) u[i].randb(nbbits);
   for (i=0; i<d; i++) v[i].randb(nbbits); 
   
  w.resize(d);
  z.resize(d);
  
  for (i=0; i<d; i++) w[i].randb(nbbits);
  for (i=0; i<d; i++) z[i].randb(nbbits);
  
  ptime.start();
  
#pragma omp parallel sections num_threads(2)
   {

#pragma omp section
     {
       for (i=0; i<d; i++)
	 v[i].addmul(u[i],alpha); 
     }

#pragma omp section
     {
       for (i=0; i<d; i++)
	 z[i].addmul(w[i],beta);
     }
   } 
   ptime.stop();

   cout << endl <<  endl << "-------- Par --------" << endl << endl;

   cout << endl << ptime << endl;
   
  
   cout << endl; 

   cout << endl <<  endl << "-------- Seq/Par --------" << endl << endl;
  
   cout << "Ratio seq/par: " << time.time()/ptime.time() << endl << endl; 
  
// #pragma omp parallel 
//   printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());



  
}


