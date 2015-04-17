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

using namespace hplll;


/***********************************************************************************

      
      Actually calls LLL with elementary lifts on ZT long et FT double 

      TODO: + Check assigment from double to Z_NR<double>
            + Use of long doubles or dpe 

      TODO: m = 1 for the moment 

      L (1 x d) and A_in d x d 

  **************************************************************************************/ 

    
  template<class ZT, class FT> int  
  lift_f_z(ZZ_mat<ZT>& U, ZZ_mat<mpz_t> L_in, ZZ_mat<FT> A_in_f, int& new_def, int def,  int target_def,
		 long shift, long increment,  int lllmethod, double delta) {

    int m,d;
    int i,j;
  
    m=L_in.getRows();
    d=L_in.getCols();

    // Toujours mpz_t
    ZZ_mat<mpz_t> L;
    L.resize(m,d);
    
    for (j=0; j<d; j++)
      L(0,j)=L_in(0,j);

    ZZ_mat<FT> Af;
    Af.resize(m+d,d);

    for (i=0; i<d; i++)
      for (j=0; j<d; j++)
	Af(m+i,j)=A_in_f(i,j);

    ZZ_mat<FT> AfT;
    AfT.resize(d,m+d);
    
    // Transform for the intermediary lifting steps 
    ZZ_mat<ZT> V;
    V.resize(d,d);
    
    ZZ_mat<FT> Vf;
    Vf.resize(d,d);

    ZZ_mat<FT> VfT;
    VfT.resize(d,d);

    // For the quit signal

    long size_of_U, size_of_V;

    // -----------
    
    Lattice<FT, FT,  matrix<Z_NR<FT> >, matrix<FP_NR<FT> > > B(Af,TRANSFORM,DEF_REDUCTION);

    // Loop
    // update def
    // get base 
    // Lift L and trucate and update Af 
    // put 

    Z_NR<mpz_t> tz;
    
    FP_NR<FT> tf;
      
    int S;

    new_def = def;
   
    
    for (S=0; S<shift; S+= increment) { 

      new_def += increment; // incrément du défaut 
      
      // Lift and truncate
     
      for (i=0; i<m; i++) 
	for (j=0; j<d; j++) {
	  tz.mul_2si(L(i,j),new_def);
	 
	  Af(i,j).getData()=tz.get_d();  // long double ?
	  
	}

     
      if (lllmethod == HLLL) {
	
	B.assign(Af);

	B.hlll(delta);
	//cout << "nblov/d: " << ((int) (((double) B.nblov)/((double) d))) << endl; 
	Af = B.getbase(); // The first row will change 
	
	Vf = B.getU();
      }
      else if (lllmethod == FPLLL) {


      	transpose(AfT,Af);
	
      	setId(VfT);

      	lllReduction(AfT, VfT, delta, 0.51, LM_FAST,FT_DEFAULT,0);

      	transpose(Af,AfT);
	
      	transpose(Vf,VfT);
	
      } 
      
      for (i=0; i<d; i++) 
	for (j=0; j<d; j++) {
	  tf = Vf(i,j).getData();  // Pour long double ou autre, vérifier et passer par set_z ? 
	  V(i,j).set_f(tf);

	}

      size_of_U = maxbitsize(U,0,d,d);
      size_of_V = maxbitsize(V,0,d,d);

      // Heuristic to check 
      if ((size_of_U + size_of_V) > 50) {

      	cerr << "**** Anomaly with the bit size of the transform (long) > 50, maybe check the value of the shift" << endl;
      	return 0;

      }
            
      matprod_in(U,V); 

      //cout << "new U bits: " << size_of_U << endl << endl;
      
      matprod_in_si(L,V);

            
    } // End main shift loop 
    
      
    return 0;

  
  };




/***********************************************************************************

     
      TODO: + Check assigment from double to Z_NR<double>
            + Use of long doubles or dpe 

  **************************************************************************************/ 
  
  template<class ZT, class FT> int lift_z(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A,  int alpha,
						long confidence_gap, long shift, long increment,
						int lllmethod, double delta) { 

    int m,d;
    int i,j;
  
    m=A.getRows();
    d=A.getCols();
    
    ZZ_mat<mpz_t> A_in;
    A_in.resize(m+d,d);
   
    // **** m=1 for the moment
    
    for (j=0; j<d; j++)
      A_in(0,j)=A(0,j);
    
    for (i=0; i<d; i++)
      A_in(m+i,i)=1;

    ZZ_mat<mpz_t> L;
    L.resize(m,d);
   
    int bitsize = maxbitsize(A,0,m,d);
  
    // For assigning the truncated basis at each step

    ZZ_mat<FT> Tf;
    Tf.resize(d,d);

    ZZ_mat<ZT> U,UT;
    U.resize(d,d);
    UT.resize(d,d);
 
    int def = -bitsize;

    int target_def = -bitsize + alpha;

    int new_def;
    
    int found=0;

   
    int intern_shift = shift;
    
    // Main loop on the shifts
    // -----------------------
    while (def < target_def) {

      HPLLL_INFO("Current default: ",def); 
      
      if ((target_def - def) <= shift) 
	intern_shift = target_def - def; 
     
      for (i=0; i<m; i++) 
     	for (j=0; j<d; j++) 
	  L(i,j)=A_in(i,j);
    
      for (i=0; i<d; i++)
     	for (j=0; j<d ; j++) 
	  Tf(i,j).getData()=A_in(m+i,j).get_d(); // long double ?

      setId(U);
      
      found=lift_f_z<ZT, FT>(U, L, Tf, new_def, def, target_def, 
			     intern_shift, increment, lllmethod, delta);

      
      
      def=new_def;

      matprod_in_si(A_in,U);

    }

    
    C=A_in;
    
    return found; // 0 here 

    
    
  }




   


/* ***********************************************

          MAIN   

   ********************************************** */



int main(int argc, char *argv[])  {
  
  
  ZZ_mat<mpz_t> A0,A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------

  int n,d;
  double delta;

  command_line_basis(A0, n, d, delta, argc, argv); 

  A.resize(n,d);
  AT.resize(d,n);
  transpose(AT,A);


  ZZ_mat<mpz_t> AR,C;
  AR.resize(1,d);
  C.resize(n,d);

  for (int j=0; j< d; j++)
    AR(0,j)=A0(0,j);

  //print2maple(AR,1,d);
  
    int start,startsec;

    Timer time;

   
    cout << "--------------  Relation type" << endl << endl; 
    start=utime();
    startsec=utimesec();
  
     
    time.start();
    lift_z<long, double>(C, AR,  maxbitsize(AR), 60, 60, 10, FPLLL, delta);
			  // long confidence_gap = 60, long shift = 200, long increment = 20,
			  // int lllmethod = FPLLL, double delta = 0.99);
    time.stop();

    start=utime()-start;
    startsec=utimesec()-startsec;
  
    
    cout << "   dimension = " << d  << endl;
    cout << "   time A: " << start/1000 << " ms" << endl;
    time.print(cout);
    cout << endl; 

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T1(C,NO_TRANSFORM,DEF_REDUCTION);
    T1.isreduced(delta-0.1);

    cout << "--------------  FPLLL WRAPPER" << endl << endl; 
    transpose(AT,A0);

    start=utime();
    startsec=utimesec();
    time.start();
    lllReduction(AT, delta, 0.501, LM_WRAPPER);
    time.stop();
    start=utime()-start;
    startsec=utimesec()-startsec;
  
    
    cout << "   dimension = " << d  << endl;
    cout << "   time B: " << start/1000 << " ms" << endl;
    time.print(cout);
 
    transpose(A,AT);
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(A,NO_TRANSFORM,DEF_REDUCTION);
    T2.isreduced(delta-0.1);

   

  return 0;
}
