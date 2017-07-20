/* LIFT LLL  

Created Jeu 16 avr 2015 14:11:32 CEST 
        

Copyright (C) 2015     Gilles Villard 

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


#ifndef HPLLL_LIFT_CC
#define HPLLL_LIFT_CC

namespace hplll {


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

    FP_NR<mpfr_t> new_quot;
    new_quot = 1.0;

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
	  Tf(i,j).get_data()=A_in(m+i,j).get_d(); // long double ?

      setId(U);
      
      found=lift_f_z<ZT, FT>(U, L, Tf, new_def, def, target_def, new_quot,
				    confidence_gap, intern_shift, increment, lllmethod, delta);

      
      
      def=new_def;

      matprod_in_si(A_in,U);

      if (found == 1) {
       
	C.resize(d,1);
	for (j=0; j<d; j++)
	  C(j,0)=A_in(m+j,0);
	
	return 1;
      }
      
    }

    // Relation bound
    // --------------

    unsigned oldprec;
    oldprec=mpfr_get_default_prec();

    mpfr_set_default_prec(2*d);

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A_in,NO_TRANSFORM,DEF_REDUCTION);

    B.householder();

    matrix<FP_NR<mpfr_t> > R;

    R=B.getR();

    FP_NR<mpfr_t> minr, tr;

    minr.abs(R(0,0));
    for (i=1; i<d; i++) {
      tr.abs(R(i,i));
      if (minr.cmp(tr) > 0) minr=tr;
    }

    cerr << endl << "** No relation found with min Rii = " << minr << endl; 

    mpfr_set_default_prec(oldprec);

    // print2maple(A_in,d+1,d);
   
    return found; // 0 here 
    
    
  }




   /***********************************************************************************

      
      Actually calls LLL with elementary lifts on ZT long et FT double 

      TODO: + Check assigment from double to Z_NR<double>
            + Use of long doubles or dpe 

      TODO: m = 1 for the moment 

      L (1 x d) and A_in d x d 

  **************************************************************************************/ 

    
  template<class ZT, class FT> int  
  lift_f_z(ZZ_mat<ZT>& U, ZZ_mat<mpz_t> L_in, ZZ_mat<FT> A_in_f, int& new_def, int def,  int target_def,
		FP_NR<mpfr_t>& new_quot,
		  long confidence_gap, long shift, long increment,
		  int lllmethod, double delta) {

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
    
    // For the end test
    // ----------------
    int found=0;

    FP_NR<mpfr_t> quot;
    

    FP_NR<mpfr_t> gap;
    gap=1.0;

    FP_NR<mpfr_t> confidence;
    // For testing 1/gap < confidence
    confidence = 1.0;
    // relié, plus petit,  au shift sur S (ex 80) 
    confidence.mul_2si(confidence,-confidence_gap); // > que increment !!! (fct taille de U ?)  

    FP_NR<mpfr_t> epsilon;
    epsilon = 10.0; // Relation to d 
    

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
	 
	  Af(i,j).get_data()=tz.get_d();  // long double ?
	  
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
	  tf = Vf(i,j).get_data();  // Pour long double ou autre, vérifier et passer par set_z ? 
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

      

      // Test
      // ----
      
      quot = new_quot;

      // 0 can be the epsilon or an artefact of integer basis reduction 
      Z_NR<mpz_t> xz;
      if (L(0,0).sgn() ==0) { // For making the gap pertinent even if 0 
	xz=1; 
      }
      else
	xz.abs(L(0,0)); 
      new_quot.set_z(xz);

      
      Z_NR<FT> tmpz,maxcol;
      
      maxcol.abs(Af(0,0));
      
      for (i=0; i<d; i++) {
       	tmpz.abs(Af(i,0));
	if (tmpz.cmp(maxcol) == 1) maxcol = tmpz;
      } 

       FP_NR<mpfr_t> xf;
       xf = maxcol.get_data(); // Double vers mpfr voir long double 
       new_quot.div(new_quot,xf);    

      gap.div(new_quot,quot);
      gap.abs(gap); 

           
      // cout << "     gap : " << gap << endl; 
      // cout << "     quot : " << new_quot << endl; 
      // cout << "     maxcol : " << maxcol << endl;
      // cout << "L: " << L(0,0) << endl;
      // cout << "Af: " << Af(0,0) << endl << endl;

      
      // Mettre avant possible
      // if (L(0,0).sgn() ==0) {
      // 	new_def = target_def;
       
      // 	return 0; 
      // }
 
      if ((gap.cmp(confidence) == -1) && (new_quot.cmp(epsilon) == -1)) {
       
	HPLLL_INFO("Candidate relation found with confidence: ",gap); 
       	return 1;
	
      } 

      
      
    } // End main shift loop 
    
      
    return found; 

  
  };



  
} // end namespace hplll


#endif 
