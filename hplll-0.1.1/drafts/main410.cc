/* Integer matrix nullspace test file  

Created Dim  7 avr 2013 16:54:03 CEST
Copyright (C) 2013      Gilles Villard 

This file is part of the hplll Library 

The hplll Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The hplll Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the hplll Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */


#include "hlll.h"
#include "matgen.h"
#include "ratio.h"
//#include "../src/nr_Z_l128.h"

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {
  
  //typedef   __int128_t  ZT;
  
  typedef   long  ZT;
  
  ZZ_mat<ZT> A; // For hpLLL 
  ZZ_mat<ZT> AT;  // fpLLL  



  int d=410;
  
  A.resize(d,d);
  
  filebuf fb;
  iostream os(&fb);

  fb.open ("pb410",ios::in);
  
  os >> A ;

  fb.close();


  double delta=0.99;

  int start,startsec;

  Timer time;

  int status;

  double t,u,v,w;
  
  Lattice<ZT, double, matrix<Z_NR<ZT> >, matrix<FP_NR<double> > > B(A,NO_TRANSFORM,SEYSEN_REDUCTION);
  
  time.start();
  status=B.hlll(delta);
  time.stop();

  start=utime()-start;
  startsec=utimesec()-startsec;
    
  cout << "   dimension = " << d  << endl;
  cout << "   time A: " << start/1000 << " ms" << endl;
  time.print(cout);

  
  //ratio<ZT>(B.getbase(),t,u,v,w);

  int dd=300;

  ZZ_mat<ZT> TB; // For hpLLL
  TB.resize(d,d);
  TB=B.getbase();
  
  ZZ_mat<ZT> T; // For hpLLL
  T.resize(d,dd);

  for (int i=0; i<d; i++)
    for (int j=0; j<dd; j++)
      T(i,j)=TB(i,j);

  ratio<ZT>(T,t,u,v,w);
  
  cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
  cout << ".. Average diagonal ratio: " << u << endl;
  cout << ".. Max diagonal ratio: " << v << endl;
  cout << ".. First vector quality: " << w << endl;

  cout << endl << "Nb swaps: " << B.nbswaps << endl;
  
  cout << endl; 
  cout << endl;
 
  //cout <<  B.getbase() << endl;

  cout << endl; 
  cout << endl;
 
  

    
  if (status ==0) {
    Lattice<ZT, mpfr_t, matrix<Z_NR<ZT> >, matrix<FP_NR<mpfr_t> > > T1(B.getbase(),NO_TRANSFORM,NO_LONG);
    T1.isreduced(delta-0.1);
  }

      
  
 
  return 0;
}
