


#include "hlll.h"

using namespace hplll;


/* --------------------------------------------------*/

int main(int argc, char *argv[])  {
 

  int m,n,i,j; 
 

  /* ----------------------------------------- */
  typedef mpz_t integer_t;

  ZZ_mat<integer_t> AT; 

  n=30;
  m=4;

  AT.resize(n,m);
  AT.gen_uniform(10);
 
  matrix<Z_NR<mpz_t> > A;
  A.resize(n,m);

  for (i=0; i<n; i++)
    for (j=0; j<m-1; j++)
      A.set(i,j,AT(i,j));

  Z_NR<mpz_t> zz;
  zz=8;

  for (i=0; i<n; i++) {
      A.set(i,m-1,AT(i,m-1));
      A(i,m-1).mul(A(i,m-1),zz);
  }


  MatrixPE<double,dpe_t> B;

  B.resize(n,m);

  for (j=0; j<m; j++)
    B.setcol(j,A.getcol(j),0,n);



  vector<FP_NR<double> >  w;
  w.resize(n);
  for (i=0; i<n; i++)
    set_z(w[i],A(i,m-1)); 

  vector<FP_NR<double> >  v;
  v.resize(n);
  for (i=0; i<n; i++)
    set_z(v[i],A(i,2)); 

  vector<double>  wf;
  wf.resize(n);
  for (i=0; i<n; i++)
    wf[i]=w[i].getData(); 

  vector<double>  vf;
  vf.resize(n);
  for (i=0; i<n; i++)
    vf[i]=v[i].getData(); 

int K=0;
 int N=16000;

 int start;

 FP_NR<double> ftmp;

  start=utime();

  for (K=0; K<N; K++)   {

    for (i=0; i<n; i++)
      //wf[i]+=ldexp(vf[i],-3);
      wf[i]+=vf[i]*8.0;
  }

 start=utime()-start;
 cout << "--- vect time:" << start << " ms" << endl;


 start=utime();
 
 for (K=0; K<N; K++)   {
   
   B.addcol(m-1,2,n);
   
 }

 start=utime()-start;
 cout << "--- mat time:" << start << " ms" << endl;

 
  //  int start,startsec,tps=0,tps_sec=0;

  /*  
  ZZ_mat<mpz_t> B;

  start = utime(); 
  startsec = utimesec(); 
  
  l1(B, A, HLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  

  cout << endl; 

  if (output) print2maple(B,n,d);
 

  cout << endl; 
  cout << "--- cputime l1 hlll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime l1 hlll is " << tps_sec << " s" << endl;
  */
  /*
  start = utime(); 
  startsec = utimesec(); 
  
  l1(B, A, FPLLL); 

  tps=utime()-start;
  tps_sec=utimesec()-startsec;

  cout << endl; 
  cout << "--- cputime l1 fplll is " << tps/1000 << " ms" << endl;
  cout << "--- cputime l1 fplll is " << tps_sec << " s" << endl;
  */




  return 0;
}
