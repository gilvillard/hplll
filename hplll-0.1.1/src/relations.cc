
/* LLL and HJLS relation algorithms 

Created Jeu  7 mar 2013 15:01:41 CET
Copyright (C) 2013,2015      Gilles Villard 

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


#include "ratio.h" 

#ifndef HPLLL_RELATIONS_CC
#define HPLLL_RELATIONS_CC

namespace hplll {


  
  /***********************************************************************************

      Relation_f  
      Restricted to doubles for the moment 
      Calls LLL with elementary lifts on ZT long et FT double 

      alpha correct bits

  **************************************************************************************/ 

  template<class ZT, class FT> 
  int relation_f(ZZ_mat<mpz_t>& C, const matrix<FP_NR<mpfr_t> > A, long alpha,
		 long confidence_gap, long shift, long increment, int lllmethod, double delta) {

    int n=A.getCols();

    ZZ_mat<mpz_t> L;
    L.resize(1,n);

    FP_NR<mpfr_t> t;
  
    for (int j=0; j<n; j++) {
      t.mul_2si( A(0,j), alpha);
      L(0,j).set_f(t);
    }

    int found;

    found=relation_f_z<ZT, FT> (C, L, alpha, confidence_gap, shift, increment, lllmethod, delta);
  
    return found;
    
  } 


  /***********************************************************************************

      Companion to relation_f  
      Relation from an integer matrix 
      Restricted to doubles for the moment 
      Calls LLL with elementary lifts on ZT long et FT double 

      alpha correct bits

      TODO: + Check assigment from double to Z_NR<double>
            + Use of long doubles or dpe 

  **************************************************************************************/ 
  
  template<class ZT, class FT> int relation_f_z(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A,  int alpha,
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
	  Tf(i,j).getData()=A_in(m+i,j).get_d(); // long double ?

      setId(U);
      
      found=detect_lift_f_z<ZT, FT>(U, L, Tf, new_def, def, target_def, new_quot,
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

      Companion to relation_f_z
      Restricted to doubles for the moment 

      Actually calls LLL with elementary lifts on ZT long et FT double 

      TODO: + Check assigment from double to Z_NR<double>
            + Use of long doubles or dpe 

      TODO: m = 1 for the moment 

      L (1 x d) and A_in d x d 

  **************************************************************************************/ 

    
  template<class ZT, class FT> int  
  detect_lift_f_z(ZZ_mat<ZT>& U, ZZ_mat<mpz_t> L_in, ZZ_mat<FT> A_in_f, int& new_def, int def,  int target_def,
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
	 
	  Af(i,j).getData()=tz.get_d();  // long double ?
	  
	}

     
      if (lllmethod == HLLL) {
	 
	B.assign(Af);

	B.hlll(delta);
	//cout << "nblov/d: " << ((int) (((double) B.nblov)/((double) d))) << endl; 
	Af = B.getbase(); // The first row will change 
	
	Vf = B.getU(); // cf conversion vers entiers implicite ??? de ZR double vers FP NR double ? bug ? 
      }
      else if (lllmethod == FPLLL) {

	{
	  // ICI 
	  ZZ_mat<mpz_t> B; 
	  B.resize(m+d,d);

	  for (int i = 0; i<m+d; i++)
	    for (int j = 0; j<d; j++)
	      B(i,j)=((long) Af(i,j).get_si());

	  double t,u,v,w;
	  ratio(B,t,u,v,w);

	  cout << endl << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	  cout << ".. Average diagonal ratio: " << u << endl;
	  cout << ".. Max diagonal ratio: " << v << endl;
	  cout << ".. First vector quality: " << w << endl;

	}


      	transpose(AfT,Af);
	
      	setId(VfT);

      	lllReduction(AfT, VfT, delta, 0.51, LM_FAST,FT_DEFAULT,0);

      	transpose(Af,AfT);
	
{
	  // ICI 
	  ZZ_mat<mpz_t> B; 
	  B.resize(m+d,d);

	  for (int i = 0; i<m+d; i++)
	    for (int j = 0; j<d; j++)
	      B(i,j)=((long) Af(i,j).get_si());

	  double t,u,v,w;
	  ratio(B,t,u,v,w);

	  cout << endl << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	  cout << ".. Average diagonal ratio: " << u << endl;
	  cout << ".. Max diagonal ratio: " << v << endl;
	  cout << ".. First vector quality: " << w << endl;

	}

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
       xf = maxcol.getData(); // Double vers mpfr voir long double 
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


  /***********************************************************************************

      Relation_lll
     
      Calls LLL with successive lifts on mpz_t and FT bases  

      alpha correct bits

   TODO: tune the parameters according to input data types 

  **************************************************************************************/ 

    
  template<class FT, class MatrixFT> int 
  relation_lll(ZZ_mat<mpz_t>& C, const matrix<FP_NR<mpfr_t> > A, long alpha,
	       long confidence_gap, long shift, int lllmethod, double delta) {

    int n=A.getCols();

    ZZ_mat<mpz_t> L;
    L.resize(1,n);

    FP_NR<mpfr_t> t;
  
    for (int j=0; j<n; j++) {
      t.mul_2si( A(0,j), alpha);
      L(0,j).set_f(t);
    }

    int found;

    found=relation_lll_z<FT, MatrixFT> (C, L, alpha, confidence_gap, shift, lllmethod, delta);

    return found;
    
  } 

  /***********************************************************************************

      Companion to relation_lll 
      Relation from an integer matrix 
      
      Calls LLL with successive lifts on mpz_t and FT bases

      alpha correct bits

      TODO: 

  **************************************************************************************/ 
  
  
  template<class FT, class MatrixFT> int  
  relation_lll_z(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A, long alpha,
		 long confidence_gap, long shift, int lllmethod, double delta) {

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

    
    int bitsize = maxbitsize(A,0,m,d);
  
    // For assigning the truncated basis at each step

    
    ZZ_mat<mpz_t> T,TT;
    
    T.resize(m+d,d);
    TT.resize(d,m+d);
 
    
    ZZ_mat<mpz_t> U,UT;
    
    U.resize(d,d);
    UT.resize(d,d);
 
    int def = -bitsize;

    int target_def = -bitsize + alpha;

    int found=0;

    FP_NR<mpfr_t> quot,new_quot,tf;
    new_quot = 1.0;  // new_quot should be bigger after the first iteration of the loop 

    FP_NR<mpfr_t> gap;
    gap=1.0;

    FP_NR<mpfr_t> confidence;
    // For testing 1/gap < confidence
    confidence = 1.0;
    // relié, plus petit,  au shift sur S (ex 80) 
    confidence.mul_2si(confidence,-confidence_gap);  // > shift  !!!  

    FP_NR<mpfr_t> epsilon;
    epsilon = 10.0; 
    
    
    Z_NR<mpz_t> tz,maxcol;

   
    Lattice<mpz_t, FT, matrix<Z_NR<mpz_t> >, MatrixFT> Bp(T,TRANSFORM,DEF_REDUCTION);
   
       
    // Main loop on the shifts
    // -----------------------
    while (def < target_def) {

      HPLLL_INFO("Current default: ",def);
       
      if ((target_def - def) <= shift) 
	def=target_def;
      else def+=shift;

      lift_truncate(T, A_in, def, shift+2*d);

      if (lllmethod == HLLL) {

	Bp.assign(T);
	
	Bp.hlll(delta);

	matprod_in(A_in,Bp.getU());
	//avec long: matprod_in_si(A_in,U);
      }

      else if (lllmethod == FPLLL) {

	transpose(TT,T);

	setId(UT);
	  
	lllReduction(TT, UT, delta, 0.51, LM_FAST,FT_DEFAULT,0);
	
	transpose(U,UT);
	
	matprod_in(A_in,U);
	//avec long: matprod_in_si(A_in,U);
	//cout << "****** sizeof U: " << maxbitsize(U,0,d,d) << endl;
      } 

      // Test
      // ----

      quot = new_quot;

      if (A_in(0,0).sgn() ==0) { // For making the gap pertinent even if 0 
	tz=1; 
      }
      else
	tz.abs(A_in(0,0)); 
      new_quot.set_z(tz);

            
      maxcol.abs(A_in(0,0));
      maxcol.mul_2si(maxcol,def);

      // cout << "L: " << A_in(0,0) << endl;
      // cout << "Af: " << maxcol << endl;
      
      for (i=0; i<d; i++) {
	tz.abs(A_in(m+i,0));
	if (tz.cmp(maxcol) ==1) maxcol = tz;
      }

      tf.set_z(maxcol);
      new_quot.div(new_quot,tf);

            
      gap.div(new_quot,quot);
      gap.abs(gap); 

            
      if ((gap.cmp(confidence) == -1) && (new_quot.cmp(epsilon) == -1)) {
       
	C.resize(d,1);
	for (j=0; j<d; j++)
	  C(j,0)=A_in(m+j,0);

	HPLLL_INFO("Candidate relation found with confidence: ",gap);
	
	return 1;
    
      } 
    
      
    } // End while 


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

   
    return found; // 0 here 

        
  };



  
} // end namespace hplll


#endif 
