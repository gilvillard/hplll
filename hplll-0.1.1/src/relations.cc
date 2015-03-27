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

    int start=utime();
    found=relation_f_z<ZT, FT> (C, L, alpha, confidence_gap, shift, increment, lllmethod, delta);
    
    start=utime()-start;

 
    cout << "     time internal: " << start/1000 << " ms" << endl;
  
    return found;
    
  } 


  /***********************************************************************************

      Companion to relation_f  
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

    cout << "     ";
    
    // Main loop on the shifts
    // -----------------------
    while (def < target_def) {

      cout << "." << flush;
      
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

    cout << endl << "** No relation found with min Rii = " << minr << endl; 

    mpfr_set_default_prec(oldprec);

   
    return found; // 0 here 
    
    
  } 
   //*************
  // Tout en double 
  // Quelques étapes en flottant à l'intérieur
  // En exact pour 
  // Pour un restart
  // FT type interne pour la matrice de passage, retournée en ZT = double  
  // avec mise à jour de def comme detect_lift 
  // **** m=1 for the moment
  // Dpe pour FT aussi, long double
  // et get_ld
  // A_in est d x d !!!!!
  // Verifier l'affectation de double à Z_NR<double>


  
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
   
    
    for (S=0; S<shift; S+= increment) {  // Limiter en borne de U  // while comme detect lift de hplll 

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

      if ((size_of_U + size_of_V) > 50) {

	cout << "**** Anomaly with the bit size of the transform (long), maybe check the value of the shift" << endl;
	return 0;

      }
            
      matprod_in(U,V); 

      //cout << "new U bits: " << size_of_U << endl << endl;
      
      matprod_in_si(L,V);

      

      // Test
      // ----
      
      quot = new_quot;
      
      Z_NR<mpz_t> xz;
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

           
      //cout << "     gap : " << gap << endl; 
      //cout << "     quot : " << new_quot << endl; 
      // cout << "     maxcol : " << maxcol << endl;
      // cout << "L: " << L(0,0) << endl;
      // cout << "Af: " << Af(0,0) << endl;

      
      // Mettre avant possible
      if (L(0,0).sgn() ==0) {
	new_def = target_def;
	return 0; 
      }
 
      if ((gap.cmp(confidence) == -1) && (new_quot.cmp(epsilon) == -1)) {
       
       	cout << endl << "     Candidate relation found with confidence " << gap << endl;  
       	return 1;
	
      } 

      
      
    } // End main shift loop 
    
      
    return found; 

  
  };


  
/*************************************************************

   Methode et limiter avec long 
   m= 1 for the moment 
   Faire une version sans troncature 
 *************************************************************/

  // Régler les paramètres en fonction des types de données
   // alpha correct bits
  
  template<class ZT, class FT, class MatrixZT, class MatrixFT> int 
  relation_lll(ZZ_mat<ZT>& C, const matrix<FP_NR<mpfr_t> > A, long alpha, long shift, int lllmethod) {

    int n=A.getCols();

    ZZ_mat<ZT> L;
    L.resize(1,n);

    FP_NR<mpfr_t> t;
  
    for (int j=0; j<n; j++) {
      t.mul_2si( A(0,j), alpha);
      L(0,j).set_f(t);
    }

    int found;

    int start=utime();
    found=relation_lll_z<ZT, FT, MatrixZT, MatrixFT> (C, L, alpha, shift, 0.99, lllmethod);
  
    start=utime()-start;
  
    cout << "   time internal: " << start/1000 << " ms" << endl;
  
    return found;
    
  } 

  //*************
  // dpe ?

  
  template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
  relation_lll_z(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int alpha, int shift, double delta, int lllmethod) {

    int m,d;
    int i,j;
  
    m=A.getRows();
    d=A.getCols();

    ZZ_mat<ZT> A_in;
    A_in.resize(m+d,d);
   
 
    // **** m=1 for the moment
    
    for (j=0; j<d; j++)
	A_in(0,j)=A(0,j);
    
    for (i=0; i<d; i++)
      A_in(m+i,i)=1;

    
    int bitsize = maxbitsize(A,0,m,d);
  
    // For assigning the truncated basis at each step

    // ICI long 
    ZZ_mat<ZT> T,TT;
    //ZZ_mat<long> T,TT;
    T.resize(m+d,d);
    TT.resize(d,m+d);
 
    // ICI long 
    ZZ_mat<ZT> U,UT;
    //ZZ_mat<long> U,UT;
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
    confidence.mul_2si(confidence,-24); // En fonction de taille de U et de dec ??? 

    FP_NR<mpfr_t> epsilon;
    epsilon = 10.0; // Relation to d 
    
    
    Z_NR<ZT> tz,maxcol;

    // ICI long 
    Lattice<ZT, FT, MatrixZT, MatrixFT> Bp(T,TRANSFORM,DEF_REDUCTION,1);
    //Lattice<long, FT, matrix<Z_NR<long> >, MatrixFT> Bp(T,TRANSFORM,DEF_REDUCTION,1);
    // Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Bp(T,TRANSFORM,DEF_REDUCTION,1);

    // Faire le produit de U ou la laisser mettre à jour ???

       
    // Main loop on the shifts
    // -----------------------
    while (def < target_def) {

      if ((target_def - def) <= shift) 
	def=target_def;
      else def+=shift;

      lift_truncate(T, A_in, def, shift+2*d);

      //cout << "****** sizeof T: " << maxbitsize(T,0,d+1,d) << "    " << sizeof(__float128) << endl;
      
      if (lllmethod == HLLL) {

	Bp.assign(T);
	
	Bp.hlll(delta);

	// ICI long 
       matprod_in(A_in,Bp.getU());
       //	matprod_in_si(A_in,Bp.getU());
      }

      else if (lllmethod == FPLLL) {

	transpose(TT,T);

	setId(UT);
	  
	lllReduction(TT, UT, delta, 0.51, LM_FAST,FT_DEFAULT,0);
	
	transpose(U,UT);
	
	// ICI long 
	matprod_in(A_in,U);
	//matprod_in_si(A_in,U);
	//cout << "****** sizeof U: " << maxbitsize(U,0,d,d) << endl;
      } 

      // Test
      // ----

      quot = new_quot;
       
      tz.abs(A_in(0,0));
      new_quot.set_z(tz);

      
      maxcol.abs(A_in(0,0));
      maxcol.mul_2si(maxcol,def);

      cout << "L: " << A_in(0,0) << endl;
      cout << "Af: " << maxcol << endl;
      
      for (i=0; i<d; i++) {
	tz.abs(A_in(m+i,0));
	if (tz.cmp(maxcol) ==1) maxcol = tz;
      }

      tf.set_z(maxcol);
      new_quot.div(new_quot,tf);

            
      gap.div(new_quot,quot);
      gap.abs(gap); 

      //cout << endl << "**  U bits: " << maxbitsize(U,0,d,d) << endl;
      cout << "     gap : " << gap << endl; 
      cout << "     quot : " << new_quot << endl; 
      cout << "     maxcol : " << maxcol << endl;

      // Mettre avant possible // Faire relation bound 
      if (A_in(0,0).sgn() ==0) {
	def = target_def;
	return 0; 
      }
      
      if ((gap.cmp(confidence) == -1) && (new_quot.cmp(epsilon) == -1)) {
       
	C.resize(1,d);
	for (j=0; j<d; j++)
	  C(0,j)=A_in(m+j,0);


	print2maple(C,1,d);
	
	cout << "Candidate relation found with confidence " << gap << endl;  
	return 1;
    
      } 
    
      
    } // End while 

    // // found = 0

    // How to do ??? One QR ??? 
    // cout << "**** There might not be relations of norm less than " << rel_bound << endl; 
  
    return found;

    
  };



  
} // end namespace hplll


#endif 
