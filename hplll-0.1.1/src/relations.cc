/* HJLS and PSLQ relation algorithms 

Created Jeu  7 mar 2013 15:01:41 CET
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



#ifndef HPLLL_RELATIONS_CC
#define HPLLL_RELATIONS_CC

namespace hplll {

  /***********************************************************************************

      Restarting 
      Restricted to doubles for the moment 
      Calls LLL with elementary lifts 

  **************************************************************************************/ 
  
  int relation_lift(ZZ_mat<mpz_t>& C, const matrix<FP_NR<mpfr_t> > A, long setprec) {

    mpfr_set_default_prec(setprec);

    int n=A.getCols();

    ZZ_mat<mpz_t> L;
    L.resize(1,n);

    FP_NR<mpfr_t> t;
  
    for (int j=0; j<n; j++) {
      t.mul_2si( A(0,j), setprec);
      L(0,j).set_f(t);
    }

    int found;

    int start=utime();
    found=relation_lift_z<mpz_t, double, matrix<FP_NR<double> > > (C, L, setprec, 0.99);
  
  
    start=utime()-start;

  
    cout << "   time internal: " << start/1000 << " ms" << endl;
  
    return found;
    
  } 

  /************************************************
     Companion integer subroutine for above one 
     Restricted to doubles for the moment 
  *************************************************/    

  template<class ZT, class FT, class MatrixFT> int  
  relation_lift_z(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int alpha, double delta) { 

    int m,d;
    int i,j;
  
    m=A.getRows();
    d=A.getCols();

    ZZ_mat<ZT> A_in;
    A_in.resize(m+d,d);

    for (i=0; i<m; i++)
      for (j=0; j<d; j++)
	A_in(i,j)=A(i,j);

    for (i=0; i<d; i++)
      A_in(m+i,i)=1;

    int bitsize = maxbitsize(A,0,m,d);
  
   
    // For assigning the exact upper lattice part at each step 
    ZZ_mat<mpz_t> L;
    L.resize(m,d);
   
    // For assigning the truncated basis at each step -- dpe_t ?? 
    ZZ_mat<double> T;
    T.resize(d+1,d);
    Z_NR<ZT> t;
   
    ZZ_mat<double> Uf;
    Uf.resize(d,d);

    ZZ_mat<long int> U;
    U.resize(d,d);

    for (i=0; i<d; i++)
      U(i,i)=1; 

    int def = -bitsize;

    int target_def = -bitsize + alpha;
   
    int new_def;

    int found;

    FP_NR<FT> rel_bound;

    Lattice<double, FT,  matrix<Z_NR<double> >,  MatrixFT> Bp(T,TRANSFORM,DEF_REDUCTION,1);
     
    // Main loop on the shifts
    // -----------------------
    while (def < target_def) {

      //int start;

      for (i=0; i<m; i++) 
	for (j=0; j<d; j++) 
	  L(i,j)=A_in(i,j);
    
      //current_shift += shift;

      //start=utime();
         
      for (i=0; i<d; i++)
	for (j=0; j<d ; j++) 
	  T(m+i,j).getData()=A_in(m+i,j).get_d(); // cf pb for assigning a double to Z_NR<double> 

      // Lattice<double, dpe_t,  matrix<Z_NR<double> >, MatrixPE<double, dpe_t> > Bp(T,TRANSFORM,DEF_REDUCTION,1);
      // FAIRE UN SET !!!!!!!!!
     
      Bp.assign(T);
      Bp.assignL(L);
     
      found = Bp.detect_lift(delta,def,target_def,new_def,maxbitsize(A_in,1,d,d),rel_bound);
      
      // start=utime()-start;
      // cout << "   time A: " << start/1000 << " ms" << endl << endl;
      // cout << " Def: " << new_def << endl;
      // start=utime();
     
     
      def=new_def;
     
      Uf=Bp.getU();

      for (i=0; i<d; i++)
	for (j=0; j<d ; j++) 
	  U(i,j).getData()=((long int) Uf(i,j).getData());
     
      matprod_in_si(A_in,U);

      if (found == 1) {
       
	C.resize(1,d);
	for (j=0; j<d; j++)
	  C(0,j)=A_in(m+j,0);
	return 1;
      }
     

    }

   
    // found = 0
    cout << "**** There might not be relations of norm less than " << rel_bound << endl; 
  
    return 0;
    
    
  } 

/*************************************************************

   Methode et limiter avec long 
   m= 1 for the moment 
   Faire une version sans troncature 
 *************************************************************/

  // Régler les paramètres en fonction des types de données
  
  template<class ZT, class FT, class MatrixZT, class MatrixFT> int 
  relation_lll(ZZ_mat<ZT>& C, const matrix<FP_NR<mpfr_t> > A, long setprec, long shift, int lllmethod=HLLL) {

    mpfr_set_default_prec(setprec);

    int n=A.getCols();

    ZZ_mat<ZT> L;
    L.resize(1,n);

    FP_NR<mpfr_t> t;
  
    for (int j=0; j<n; j++) {
      t.mul_2si( A(0,j), setprec);
      L(0,j).set_f(t);
    }

    int found;

    int start=utime();
    found=relation_lll_z<ZT, FT, MatrixZT, MatrixFT> (C, L, setprec, shift, 0.99, lllmethod);
  
  
    start=utime()-start;

  
    cout << "   time internal: " << start/1000 << " ms" << endl;
  
    return found;
    
  } 

  //*************
  
  template<class ZT, class FT, class MatrixZT, class MatrixFT> int  
  relation_lll_z(ZZ_mat<ZT>& C, ZZ_mat<ZT> A, int alpha, int shift, double delta, int lllmethod=HLLL) {

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
    //ZZ_mat<ZT> T,TT;
    ZZ_mat<long> T,TT;
    T.resize(m+d,d);
    TT.resize(d,m+d);
 
    // ICI long 
    //ZZ_mat<ZT> U,UT;
    ZZ_mat<long> U,UT;
    U.resize(d,d);
    UT.resize(d,d);
 
    int def = -bitsize;

    int target_def = -bitsize + alpha;

    int found;

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
    //Lattice<ZT, FT, MatrixZT, MatrixFT> Bp(T,TRANSFORM,DEF_REDUCTION,1);
    Lattice<long, FT, matrix<Z_NR<long> >, MatrixFT> Bp(T,TRANSFORM,DEF_REDUCTION,1);
    // Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > Bp(T,TRANSFORM,DEF_REDUCTION,1);

    // Faire le produi de U ou la laisser mettre à jour ???

    
    
    // Main loop on the shifts
    // -----------------------
    while (def < target_def) {

      if ((target_def - def) < shift) def=target_def;
      else def+=shift;

      lift_truncate(T,A_in,def,shift+d);

      cout << "****** sizeof T: " << maxbitsize(T,0,d+1,d) << endl;
      
      if (lllmethod == HLLL) {
	      
	Bp.assign(T);
	
	Bp.hlll(delta);

	// ICI long 
	//matprod_in(A_in,Bp.getU());
	matprod_in_si(A_in,Bp.getU());
      }

      else if (lllmethod == FPLLL) {

	transpose(TT,T);

	setId(UT);
	
	lllReduction(TT, UT, delta, 0.51, LM_FAST,FT_DEFAULT,0);

	transpose(U,UT);
	  
	// ICI long 
	//matprod_in(A_in,U);
	matprod_in_si(A_in,U);
	cout << "****** sizeof U: " << maxbitsize(U,0,d,d) << endl;

      } 
      
      // Test
      // ----

      quot = new_quot;
       
      tz.abs(A_in(0,0));
      new_quot.set_z(tz);

      
      maxcol.abs(A_in(0,0));
      maxcol.mul_2si(maxcol,def);

      for (i=0; i<d; i++) {
	tz.abs(A_in(m+i,0));
	if (tz.cmp(maxcol) ==1) maxcol = tz;
      }

      tf.set_z(maxcol);
      new_quot.div(new_quot,tf);

            
      gap.div(new_quot,quot);
      gap.abs(gap); 

      if ((gap.cmp(confidence) == -1) && (new_quot.cmp(epsilon) == -1)) {
       
	C.resize(1,d);
	for (j=0; j<d; j++)
	  C(0,j)=A_in(m+j,0);


	print2maple(C,1,d);
	
      cout << "Candidate relation found with confidence " << gap << endl;  
      return 1;
    
      found=0;
      } 
      
    } // End while 

    // // found = 0

    // How to do ??? One QR ??? 
    // cout << "**** There might not be relations of norm less than " << rel_bound << endl; 
  
    return found;

    
  };

  
// ********   ATTENTION   ******************************
//
// Template RT but written for mpfr .....
//
/* *************************************************************

   COMPUTING INTEGER RELATIONS FOR A FLOATING POINT MATRIX 

   Method 1: HJLS, starting from the matrix d x n 
              using the identiy for calling decomp with dec > 0 

   Method 2: 1) compute an associated floating FGAS 
             2) go through the FGAS decomposition      

  
   d x n 
   Hence nullity (returned) is <= nbrel < n-d */

   

/***********************************************************************************

   HJLS 

   Starting from the matrix d x n using the identiy for calling decomp with dec > 0 

**************************************************************************************/ 


// Not templated w.r.t. other matrices for the moment at this upper level 
// Put the integer matrices 
template<class RT, class FT, class MatrixFT>  int relations_hjls(ZZ_mat<mpz_t>& C, const matrix<FP_NR<RT> >& A, long nbrel, long setprec) {

  int d,n;

  d=A.getRows();
  n=A.getCols();

  mpfr_set_default_prec(setprec); // Same as for decomp? 

  matrix<FP_NR<RT> > B;

  B.resize(n,n+d);

  int i,j;

  for (i=0; i<n; i++) {
    for (j=0; j<d; j++)
      B(i,j)=A(j,i);  // = is allowed, set is not needed? 
    B(i,d+i)=1;
  }



  Fgas<RT, mpz_t, FT,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, MatrixFT > F(B, NO_TRANSFORM, setprec,d);
     
  int flag_decomp;

 
  flag_decomp=F.decomp(1.1547005384, nbrel, LONG_MAX);
  //flag_decomp=F.decomp(1.4142135624, nbrel, LONG_MAX);

  if (flag_decomp==-1) {
    return -1;
  }


  // The nullspace C should be in the last n+d-rel columns of V=Transpose(U^(-1))
  // See analogous procedure in nullspace indirect (on other entry types) 
  // ******************************************************************************

  ZZ_mat<mpz_t> UT;
  UT.resize(n+d,n+d);
  transpose(UT,F.getV());

  C.resize(n,nbrel);

  for (i=0; i<n; i++)
    for (j=0; j<nbrel; j++) 
      C(i,j)=UT(n+d-nbrel+j,d+i); 

  // Check of the nullity 
  // Coul maybe not be used id decomp has succeeded?  
  // ***********************************************

  matrix<FP_NR<RT> > Cfp, ZM;

  Cfp.resize(n,nbrel);
  for (i=0; i<n; i++)
    for (j=0; j<nbrel; j++)
      set_z(Cfp(i,j),C(i,j));

  ZM.resize(d,nbrel);

  

  matprod(ZM,A,Cfp); // Should be epsilon-zero 
  

  /* Tuning of the epsilon = "zero" value for the test, about the floating point precision */
  /* ***************************************************************************************/

  FP_NR<RT> epsilon,tmp1;

  epsilon=F.epsilon;   // Should/could use a different epsilon? Yes if the precision of the input A is greater 


  /* Zero test    */
  /****************/

  int nullity=0;
  bool zeroedcol;

  
  for (j=0; j<nbrel; j++) {
       
    zeroedcol=true;

    for (i=0; i<d; i++) {
    
      tmp1.abs(ZM(i,j));  
      if (tmp1.cmp(epsilon) > 0)   
	zeroedcol=false;
    }
    if (zeroedcol) 
      nullity+=1;   
  }


  if (nullity != nbrel) {
    cout << " *******  ERROR: Input matrix non-full row rank or problem with the decomposition/precision" << endl; 
    cout << " *******          nullity = " << nullity << " <> " << nbrel << endl;  
    return -1;
  } 


  return(nbrel);

} // End body of FP relations HJLS  



   

/***********************************************************************************

   RESTARTING HJLS 

   Starting from the matrix d x n using the identiy for calling decomp with dec > 0 

**************************************************************************************/ 


// Not templated w.r.t. other matrices for the moment at this upper level 
// Put the integer matrices 
template<class RT, class FT, class MatrixFT>  int restarting_hjls(ZZ_mat<mpz_t>& C, const matrix<FP_NR<RT> >& A, long nbrel, long setprec) {

  int d,n;

  int nblov=0;

  d=A.getRows();
  n=A.getCols();

  int N;
  N=n+d;

  int tps=0; 
  int start; 


  FP_NR<RT> tr,acc; // Tmp  

  int i,j,k;
 
  mpfr_set_default_prec(setprec); // Same as for decomp?

  // Main real matrix 
  // ----------------
 
  // Main FGAS 
  matrix<FP_NR<RT> > B;
  B.resize(n,N);

  for (i=0; i<n; i++) {
    for (j=0; j<d; j++)
      B(i,j)=A(j,i);  // = is allowed, set is not needed? 
    B(i,d+i)=1;
  }

  matrix<FP_NR<RT> > Res;
  Res.resize(n,d);
 
  matrix<FP_NR<RT> > H;
  H.resize(n,N);

 
  // Float matrices 
  // --------------
  MatrixFT Bft;

  Bft.resize(n,N);

  MatrixFT Resft;

  Resft.resize(n,d);


  // Transform and inverse transpose 
  // -------------------------------

  // Global 
  ZZ_mat<mpz_t> U; 

  U.resize(N,N);
  for (i=0; i<N; i++)
    U(i,i)=1;

  // Local decomposition transformation 
  //ZZ_mat<mpz_t> FU; 
  ZZ_mat<long> FU;

  // Temporary FU in real format 
  matrix<FP_NR<RT> > FUr;
  FUr.resize(N,N);

  // Global inverse transpose transformation 
  ZZ_mat<mpz_t> V; 
  
  V.resize(N,N);
  for (i=0; i<N; i++)
    V(i,i)=1;
  
  // For later zero-epsilon test at high precision 
  FP_NR<RT> epsilon; 


  // Initial size reduction, especially for input cases
  // with very different magnitudes
  // And high precision FGAS, first loop preparation 

 

  Fgas<RT, mpz_t, RT,  matrix<FP_NR<RT> >, matrix<Z_NR<mpz_t> >, matrix<FP_NR<RT> > > G(B,NO_TRANSFORM,setprec);

  // ICI ce shift only for 1010 pb understanding 
  //G.shift_epsilon(220);
  epsilon=G.epsilon; 

  for (j=0; j<N; j++) 
    G.hsizereduce(j);

  matprod_in(V,G.getV());
  B=G.getfgas();
  H=G.getR();

  for (j=0; j<N; j++)
    Bft.setcol(j,H.getcol(j),0,n);
   


  // Residue update for next round 
  // -----------------------------

  for (i=0; i<n; i++)
    for (j=0; j<d; j++) {
      
      acc=0.0;
      for (k=0; k<n; k++) {
	set_z(tr,V(d+k,d+i));
	tr.mul(tr,A(j,k));  
	acc.add(acc,tr);
      }
      Res.set(i,j,acc);
    } 
  
  // Residue normalization  
  acc=0.0;
  for (i=0; i<n; i++)
    for (j=0; j<d; j++) {
      tr.abs(Res(i,j));
      if (tr.cmp(acc)>0) acc=tr;
    }
  
  for (j=0; j<d; j++)
    Res.div(j,0,Res.getcol(j,0),acc,n);
  
  for (j=0; j<d; j++)
    Resft.setcol(j,Res.getcol(j),0,n);
  

  // Main loop preparation 
  // ---------------------

  FP_NR<RT> r1;

  int flag_decomp;
  int nbloops=0;


  bool notfound=1;
  long internal_prec;
  internal_prec=(Resft.get(0,0)).getprec();


  // Low precision FGAS
  // ------------------

  Fgas<FT, long, FT,  MatrixFT, matrix<Z_NR<long> >, MatrixFT > F(Bft, Resft, TRANSFORM, internal_prec, d);

   
  F.shift_epsilon(16);
  F.shift_testU(4);
  // ICI que pour 1010 
  //F.shift_epsilon(12);
  //F.shift_testU(5);

  // ****************************************
  // Main loop on low precision decomposition
  // ****************************************


  while (notfound) {
  
    nbloops+=1;

    // FGAS decomposition  
    // ------------------

    F.set(Bft, Resft);
   
    start=utime();
 
    // ICI   restarting 
    flag_decomp=F.decomp(1.1547005384, nbrel, LONG_MAX);
    //flag_decomp=F.decomp(1.4142135624, nbrel, LONG_MAX);

    mpfr_set_default_prec(setprec); // Should be only an internal change above (but changes everything) 

    start = utime()-start;
    tps+=start;
    
    nblov+=F.nblov;

    
    
    // If relation found 
    // -----------------
    if (flag_decomp==0) {
    
      matprod_in_si(V,F.getV());
      notfound=0;

      // ICI 
      cout << endl << " WITH SMALL DECOMP ENDING " << endl << endl; 

      // ICI 
      print2maple(V,n+d,n+d);

    }
    else { // Relation not found, will continue 
           // ---------------------------------

      FU=F.getU();
      
      for (i=0; i<N; i++) 
	for (j=0; j<N; j++) {
	  set_z(tr,FU(i,j));
	  FUr.set(i,j,tr);
	}

     
      // Accumulate in B 
      // ---------------
      
      matprod_in(B,FUr);
      
      // V accumulation 
      // --------------
      
      matprod_in_si(V,F.getV());

      // and size-reduction for better numerical quality 
      // ----------------------------------------------- 

      G.set(B);

      for (j=0; j<N; j++) {
	G.hsizereduce(j);
      }

      matprod_in(V,G.getV());
      B=G.getfgas();

      H=G.getR();

      for (j=0; j<N; j++)
      Bft.setcol(j,H.getcol(j),0,n);

      
      // Residue update for next round 
      // -----------------------------

      for (i=0; i<n; i++)
	for (j=0; j<d; j++) {
	  
	  acc=0.0;
	  for (k=0; k<n; k++) {
	    set_z(tr,V(d+k,d+i));
	    tr.mul(tr,A(j,k));  
	    acc.add(acc,tr);
	  }
	  Res.set(i,j,acc);
	} 
     
      // Early zero detection 
      // --------------------
      // TODO : also put heuristic for global precision exhaustion 
      r1=0.0;
      for (j=0; j<d; j++) {
	tr.abs(Res.get(n-1,j)); // Change for several relations and ldim 
	if (tr.cmp(r1) > 0) r1=tr;
      }

      if (epsilon.cmp(r1) >0) {
       
	notfound=0; 
	// ICI 
	cout << endl << " WITH RESIDUE ENDING " << endl << endl; 
      }
 
      // Residue normalization  
      acc=0.0;
      for (i=0; i<n; i++)
	for (j=0; j<d; j++) {
	  tr.abs(Res(i,j));
	  if (tr.cmp(acc)>0) acc=tr;
	}
    
      for (j=0; j<d; j++)
	Res.div(j,0,Res.getcol(j,0),acc,n);

      for (j=0; j<d; j++)
	Resft.setcol(j,Res.getcol(j),0,n);

 
    } // End else not found, continue 

  } // End main loop on successive decompositions  

 

  // ICI 
  cout << endl << "   Part time: " << tps/1000 << " ms" << endl;
  cout << endl << "   Nb iterations " << nblov << endl;

  // The nullspace C should be in the last n+d-rel columns of V=Transpose(U^(-1))
  // See analogous procedure in nullspace indirect (on other entry types) 
  // ******************************************************************************
 
 

  ZZ_mat<mpz_t> UT;
  UT.resize(n+d,n+d);
  transpose(UT,V);
  
  C.resize(n,nbrel);
  
  for (i=0; i<n; i++)
    for (j=0; j<nbrel; j++) 
      C(i,j)=UT(n+d-nbrel+j,d+i); 
  
  
  // Check of the nullity 
  // Coul maybe not be used id decomp has succeeded?  
  // ***********************************************

  
  matrix<FP_NR<RT> > Cfp, ZM;

  Cfp.resize(n,nbrel);
  for (i=0; i<n; i++)
    for (j=0; j<nbrel; j++)
      set_z(Cfp(i,j),C(i,j));

  ZM.resize(d,nbrel);

  matprod(ZM,A,Cfp); // Should be epsilon-zero 
  
  /* Tuning of the epsilon = "zero" value for the test, about the floating point precision */
  /* ***************************************************************************************/

  
  FP_NR<RT> tmp1;

  /* Zero test    */
  /****************/
  
  int nullity=0;
  bool zeroedcol;

  for (j=0; j<nbrel; j++) {
       
    zeroedcol=true;

    for (i=0; i<d; i++) {
    
      tmp1.abs(ZM(i,j));     
      if (tmp1.cmp(epsilon) > 0)   
	zeroedcol=false;
    }
    if (zeroedcol) 
      nullity+=1;   
  }

 
  if (nullity != 1) {
    cout << " ***  ERROR: Input matrix non-full row rank or problem with the decomposition/precision" << endl; 
    cout << " ***          nullity = " << nullity << " <> " << nbrel << endl;  

    for (j=0; j<nbrel; j++) {
      
      tr=0.0;
      for (i=0; i<d; i++) {
	tmp1.abs(ZM(i,j));
	if (tmp1.cmp(tr) >0) tr=tmp1;
      }
    }
    cout << " ***          max residue = " << tr << "  >  epsilon = " << epsilon << endl;
    
    return -1;
  } 

  return(nbrel);
  
} // End body of FP relations HJLS  

} // end namespace hplll


#endif 
