/* Householder LLL 

Created Mar 18 jan 2011 18:10:25 CET 
Copyright (C) 2011, 2012, 2013      Gilles Villard 

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

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef HPLLL_HLLL_CC
#define HPLLL_HLLL_CC


// MPFR setprec : isreduced and cond are re-initializing the floating point matrices 

namespace hplll {   

  
template<class ZT,class FT, class MatrixZT, class MatrixFT>  int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::hlll(double delta, bool verbose) { 

 
  
  if (verbose) 
    verboseDepth-=1;
  
  int kappa=1,i;
  int prevkappa=-1; // For the looping test between tow indices 
  vector<FP_NR<FT> >  prevR(d);

  Timer cpu_discovered,cpu_tot;
  cpu_discovered.clear();
  cpu_tot.clear();
 

  FP_NR<FT> newt; //testaccu;

  FP_NR<FT> deltab,lovtest;
  deltab=delta;   // TO SEE 

  FP_NR<FT> tmpswap;

  FP_NR<FT> s,sn; // newt test 

  int flag_reduce=0; // No convergence in reduce 


  for (i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}

  // **********
  // int kmax=2;

  // Timer time;

  // time.start();

  // int nb=0;
  
  //************
  
  cpu_discovered.start();

  while ((kappa < kappa_max) && (nblov < nblov_max)) 
    {

     if (verbose) 
       if (verboseDepth >= 0) {
	 if (((nblov%800000)==0) && (nblov > 0))   cout << nblov << " tests" << endl; 
       }
      
      if (kappa == 1) { 
	for (i=0; i<d; i++) kappamin[i]=min(kappamin[i],0);
 
	householder_r(0);  
	householder_v(0);
	
      } 
      else for (i=kappa+1; i<d; i++) kappamin[i]=min(kappamin[i],kappa); // lowest visit after kappa 

      
      if (seysen_flag < 1) 
      	flag_reduce=hsizereduce(kappa);
      else
	flag_reduce=seysenreduce(kappa);
     
                  
      if (flag_reduce==-1) {
	if (verbose) 
	  verboseDepth+=1;  // Segment à commenter 
	return(-1);
      }

      //newt=normB2[kappa];  // The newt test was like that in old non exp 
      //for (i=0; i<=kappa-2; i++) newt.submul(R.get(i,kappa),R.get(i,kappa));

      lovtest.mul(R.get(kappa-1,kappa-1),R.get(kappa-1,kappa-1));
      lovtest.mul(deltab,lovtest);

      nblov+=1;
    
      fp_norm(s,R.getcol(kappa,kappa),structure[kappa]+1-kappa);
      
      s.mul(s,s);

      
      
      sn.mul(R.get(kappa-1,kappa),R.get(kappa-1,kappa));
      newt.add(s,sn);
 

      // ****************
      //   UP    UP    UP 
      // ****************
    
      if (lovtest <= newt) {

	
	// Cpu time, new kappa discovered
	// ------------------------------
	
	if (((kappa) < d) && (col_kept[kappa+1] == 0)) {
	  cpu_discovered.stop();

	  cpu_tot+=cpu_discovered;
	  
	  if (verbose) 
	    if (verboseDepth >= 0) {
	      if (kappa < d-1) cout << "Discovering vector " << kappa +2 << "/" << d << endl;
	      else cout << endl; 
	      cout << "     Phase-: " << cpu_discovered << endl;
	      cout << "     Total: " << cpu_tot << endl;
	    }

	 
	  cpu_discovered.start(); 
	} 
							      
	
	householder_v(kappa);   // The first part of the orthogonalization is available 

	
	// Heuristique precision check : when R(kappa-1,kappa-1) increases in a 2x2 up and down  
	// ------------------------------------------------------------------------------------
	if (prevkappa==kappa+1) {  
	  FP_NR<FT> t;
	  t.abs(R.get(kappa,kappa));
       
	  if (t > prevR[kappa]) {
	
	    cout << " **** #tests = " << nblov << " **** Anomaly: the norm increases for kappa = " << kappa << endl;

	    if (verbose) 
	      verboseDepth+=1;   // Segment à commenter 
	    return -1;
	  }

	}

	

	descendu[kappa]=0;
	if (kappa==1) descendu[0]=0;
	
	prevkappa=kappa; 
	prevR[kappa].abs(R.get(kappa,kappa)); 
	
	kappa+=1; 



      } // End up 

      // ****************
      //   DOWN   DOWN 
      // ****************

      else {

	if (kappa==1) descendu[0]=0 ;
	else descendu[kappa-1]=descendu[kappa]+1;
	descendu[kappa]=0;

	nbswaps+=1;
       
	B.colswap(kappa-1,kappa);
	
	if (transf) U.colswap(kappa-1,kappa);
	
	Bfp.colswap(kappa-1,kappa);

	structure[kappa-1]=structure[kappa];

	VR.colswap(kappa-1,kappa);
 
	tmpswap=normB2[kappa];  normB2[kappa]=normB2[kappa-1]; normB2[kappa-1]=tmpswap;

	prevkappa=kappa; 
	kappa=max(kappa-1,1);  

      }

    } // End main LLL loop 

  if (verbose) 
    verboseDepth+=1;  // Segment à commenter 
  return 0;
  
};



/* -------------------------------------------------------------------------
   Seysen size reduction 

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm

   ------------------------------------------------------------------------- */


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::seysenreduce(int kappa) { 

  
  nmaxkappa=structure[kappa]+1;

  FP_NR<FT> approx;
  approx=0.01;
   
  FP_NR<FT>  qq;
  
  FP_NR<FT> x,t;

  vector<FP_NR<FT> > vectx(kappa);  
  
  vector<FP_NR<FT> > tmpcolR(kappa);  
 
  int i,k,w=0;

  bool nonstop=1;
  bool somedone=0;

  vector<bool> bounded(kappa);
  
  int restdim=0; // Remaining dimension after the current block 

  
  householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 


  FP_NR<FT> theta,eps;
  

  int bdim,ld,tdig,indexdec;
  
  while (nonstop) {  // LOOP COLUMN CONVERGENCE

  
  
    // Jeu 12 mai 2016 12:55:13 CEST 
    theta.sqrt(normB2[kappa]);
    eps=0.00000000001; // To tune for theta depending on the precision
    theta.mul(theta,eps);
  
    w++;

    somedone = 0;

        
    ld=1; indexdec=0; // Décalage d'indice
 

    while (ld <=kappa) {
      
      tdig=(kappa/ld)%2;
      
      if (tdig==0) bdim =0; else bdim = ld;

      // -----------------------------------------------------
      // Boucle sur la partie de la colonne correspond au bloc 
      // -----------------------------------------------------

      // vectxz rounding of vectx
      //  column -  (prev col) * vectxz (xz ré-utilisé) 
      // On peut travailler sur place en remontant dans la colonne kappa de R 

      // On calcule vectx et on arrondit au fur et à mesure

      restdim=kappa-indexdec-bdim;

		       
      for (i=kappa-1-indexdec; i>=restdim; i--) 
	tmpcolR[i]=R.get(i,kappa);

      compteur =0;
      
      for (i=kappa-1-indexdec; i>=restdim; i--){
	
	vectx[i].div(tmpcolR[i],R.get(i,i)); // Faire après le test, pas besoin si bounded ?

	qq.abs(R.get(i,i)); 
	qq.add(qq,theta); // À optimiser 
	qq.div(tmpcolR[i],qq);
	qq.abs(qq);
	//qq.abs(vectx[i]);
	
	if (qq.cmp(0.501) == 1) {
	  bounded[i]=0;
	  somedone=1;
	}
	else bounded[i]=1;

		
	// Faire une opération vectorielle 
	for (k=restdim; k<i; k++) 
	  tmpcolR[k].submul(R.get(k,i),vectx[i]);
	 
      } // end calcul de la transfo 
       

      // Et on applique la transformation  
      // --------------------------------

      
      // Sequential
      // ----------

      if ((num_S==1) || (bdim < par_reduce_threshold))  
	seysen_update(vectx, kappa, kappa-1-indexdec, restdim,  bounded);

      
      // Parallel     
      // --------
      else 
	{
	  seysen_update_R(vectx, kappa, kappa-1-indexdec,  restdim,  bounded);
	  pseysen_update_B(kappa, kappa-1-indexdec,  restdim, vectx, bounded, num_S);
	}
      
      indexdec+=bdim;     
      ld=ld*2;

    } // End loop on log blocks 

    

    if (somedone) {
       
      col_kept[kappa]=0;

      t.mul(approx,normB2[kappa]);

      householder_r(kappa);  
     
    }
    else {
      
      nonstop=0;
     
     }

  } // end while 

  return somedone;

};

/* -------------------------------------------------------------------------
   Seysen sequential update of R and B   

   ------------------------------------------------------------------------- */

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline bool 
  Lattice<ZT,FT, MatrixZT, MatrixFT>::seysen_update( vector<FP_NR<FT> >& vectx, int kappa, int from_i, int restdim, vector<bool> bounded) { 


    int i;
    
    FP_NR<FT>  xf; // Old x

    long expo,lx;

    int nmax;
    
    Z_NR<ZT>  xz;
    
    for (i=from_i; i>= restdim; i--){

	vectx[i].rnd(vectx[i]);
	xf=vectx[i]; 

	
	if (bounded[i]==0) {

	  	
	  lx = xf.get_si_exp(expo);

	  nmax=structure[i]+1;
	  
	  // Cf fplll 
	  // Long case 
	  if (expo == 0) {
	    
	    if (lx == 1) {

	      R.subcol(kappa,i,restdim);
	      
	      B.subcol(kappa,i,nmax);
	      
	      if (transf) 
		U.subcol(kappa,i,min(d,nmax));
	      
	    } 
	    else if (lx == -1) {

	     
	      R.addcol(kappa,i,restdim);
	      
	      B.addcol(kappa,i,nmax);
	      
	      if (transf) 
		U.addcol(kappa,i,min(d,nmax));

	    } 
	    else { 
	      
	      if (fast_long_flag == 1) {
	      
	       	R.submulcol(kappa,i,xf,restdim);
		
	       	B.addmulcol_si(kappa,i,-lx,nmax);
	       	if (transf)  
	       	  U.addmulcol_si(kappa,i,-lx,min(d,nmax));
		
	      } // end fast_long
	      else {
		
	       	set_f(xz,x);  
	      
	       	R.submulcol(kappa,i,xf,restdim);	
	       	B.submulcol(kappa,i,xz,nmax);
	       	if (transf)  
	       	  U.submulcol(kappa,i,xz,min(d,nmax));
	      }	      
	    } 
	    
	  } // end expo == 0 
	  else {  // expo <> 0 

	   
	    if (fast_long_flag == 1) {
	    
	      R.submulcol(kappa,i,xf,restdim);
	      
	      B.addmulcol_si_2exp(kappa,i,-lx,expo,nmax);
	      if (transf)  
		U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));
	    
	    } // end fast_long
	    else {
	   
	      set_f(xz,x);  
	  
	      R.submulcol(kappa,i,xf,restdim);
	      
	      B.submulcol(kappa,i,xz,nmax);
	      if (transf)  
		U.submulcol(kappa,i,xz,min(d,nmax));
	  
	    } // end no long


	    
	  } // end expo <> 0 
	} // Non zero combination 

      } // end application de la transformation 
    

    return 0;

  }


  
/* -------------------------------------------------------------------------
   Seysen update of R  / TO DO parallel 

   ------------------------------------------------------------------------- */

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline bool 
  Lattice<ZT,FT, MatrixZT, MatrixFT>::seysen_update_R(vector<FP_NR<FT> >& vectx, int kappa, int from_i, int restdim, vector<bool> bounded) { 


    int i;
    
    FP_NR<FT>  xf; // Old x

    long expo,lx;
    
    
    for (i=from_i; i>= restdim; i--){

      vectx[i].rnd(vectx[i]);  // Attention si en parallèle 
      xf=vectx[i]; 

	
      if (bounded[i]==0) {

	  	
	lx = xf.get_si_exp(expo);

	// Cf fplll 
	// Long case 
	if (expo == 0) {
	    
	  if (lx == 1) {

	    R.subcol(kappa,i,restdim);
	      
	  } 
	  else if (lx == -1) {

	     
	    R.addcol(kappa,i,restdim);
	      
	  } 
	  else { 
	      
	    if (fast_long_flag == 1) {
	      
	      R.submulcol(kappa,i,xf,restdim);
		
	    } // end fast_long
	    else {
		
	      
	      R.submulcol(kappa,i,xf,restdim);	
	       	
	    }	      
	  } 
	    
	} // end expo == 0 
	else {  // expo <> 0 

	  R.submulcol(kappa,i,xf,restdim);

	} // end expo <> 0 
      } // Non zero combination 

    } // end application de la transformation 
    

    return 0;

  }


 

    
  /* -------------------------------------------------------------------------
   Seysen parallel update of B  

   ------------------------------------------------------------------------- */

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline bool 
  Lattice<ZT,FT, MatrixZT, MatrixFT>::pseysen_update_B(int kappa, int from_i, int restdim,
							vector<FP_NR<FT> > vectx, vector<bool> bounded, int S) { 

    
    int nmax;
      
    nmax=structure[from_i];
    
    for (int i=from_i+1; i>=restdim; i--)  
      nmax=max(nmax,structure[i]);
    
    nmax +=1;

    // Parallel loop
    // -------------
    
#ifdef _OPENMP
#pragma omp parallel for 
#endif

    for (int l=0; l<S; l++) {

	int i;
       
	FP_NR<FT>  xf; // Old x
      
	long expo,lx;
      
	Z_NR<ZT>  xz;

	int nloc;
	nloc=nmax/S;
	if ((l+1) <= nmax%S)
	  nloc+=1;
      
	int ibeg;
	ibeg=(nmax/S)*l;
	if ((l+1) <= nmax%S)
	  ibeg+=l;
	else
	  ibeg+=nmax%S;

      
	for (i=from_i; i>= restdim; i--){
    
	  xf=vectx[i]; 
      
	
	  if (bounded[i]==0) {
	  
	  
	    lx = xf.get_si_exp(expo);
	  
	    // Cf fplll 
	    // Long case 
	    if (expo == 0) {
	    
	      if (lx == 1) {

	     
		B.psubcol(kappa,i,ibeg,nloc);
	      
		// if (transf)   TODO EN PARALLELE CF LONGUEUR
		// 	U.subcol(kappa,i,min(d,nmax));
	      
	      } 
	      else if (lx == -1) {

		B.paddcol(kappa,i,ibeg,nloc);
	      
		// if (transf) 
		// 	U.addcol(kappa,i,min(d,nmax));

	      } 
	      else { 
	      
		if (fast_long_flag == 1) {
	      
	       
		  B.paddmulcol_si(kappa,i,-lx,ibeg,nloc);
		
		  // if (transf)  
		  //   U.addmulcol_si(kappa,i,-lx,min(d,nmax));
		
		} // end fast_long
		else {
		
		  set_f(xz,x);  
	      
		  B.psubmulcol(kappa,i,xz,ibeg,nloc);
		
		  // if (transf)  
		  //   U.submulcol(kappa,i,xz,min(d,nmax));
		}	      
	      } 
	    
	    } // end expo == 0 
	    else {  // expo <> 0 

	   
	      if (fast_long_flag == 1) {
	    
		B.paddmulcol_si_2exp(kappa,i,-lx,expo,ibeg,nloc);
	      
		// if (transf)  
		// 	U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));
	    
	      } // end fast_long
	      else {
	   
		set_f(xz,x);  
	  
		B.psubmulcol(kappa,i,xz,ibeg,nloc);
	      
		// if (transf)  
		// 	U.submulcol(kappa,i,xz,min(d,nmax));
	  
	      } // end no long

	      
	    
	    } // end expo <> 0 
	  } // Non zero combination 

	} // end application de la transformation 

    } // Parallel loop 

    return 0;

  }
  

 


/* -------------------------------------------------------------------------
   Size reduction 

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm

   ------------------------------------------------------------------------- */


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::hsizereduce(int kappa, int fromk) { 

  
  nmaxkappa=structure[kappa]+1;

  vector<FP_NR<FT> > vectx;
  vectx.resize(n);
  
  FP_NR<FT> approx;
  
  approx=0.1;


  FP_NR<FT> t;
 
  int i,w=0;

  bool nonstop=1;
  bool prev_nonstop = 1;
  
  bool somedone=0;

 
  
  householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous

  
      
  // While loop for the norm decrease
  // --------------------------------
  

  int startposition;
  if (fromk > 0) 
    startposition = min(kappa-1,fromk);
  else 
    startposition = kappa-1;

    somedone = 1;
  
  //while (nonstop) {
  while (somedone == 1) {

    
    w++;


    // Sequential
    // ----------

   
     
    if ((num_S==1) || ((kappa-1) < par_reduce_threshold))  
      somedone=size_update(kappa, startposition, 0);

    // Parallel
    // --------
    else {

     
      somedone=size_update_R(vectx, kappa, startposition, 0); 
      
      psize_update_B(kappa, startposition, 0, vectx, num_S);
    } 


          
    if (somedone) {
      
      col_kept[kappa]=0;
      
      t.mul(approx,normB2[kappa]);
      
      householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 

      nonstop = (normB2[kappa] < t);  // ne baisse quasiment plus ?

      // Heuristic test
      // The norm is not increasing for several steps
      // This may happen exceptionnaly in correct cases with mu_ij = 1/2 exactly
      //  (and alternates between to values 1/2 and -1/2 in floating point 
      // Or happen when not enough precision
      // Hence here: test of size reduction, if yes then exit if no then return -1 
      if ((prev_nonstop ==0) && (nonstop == 0)) {
      	
      	  FP_NR<FT> one;
      	  one = 1.0;
	  
      	  FP_NR<FT> theta,eps;
      	  //theta = 0.0000001;
      	  //theta.mul(theta,R.get(kappa,kappa));
	  // Jeu 12 mai 2016 12:55:13 CEST - R.get(kappa,kappa) may not be relevant (no hoseholder_v) 
	  theta.sqrt(normB2[kappa]);
	  eps=0.00000000001; // To tune for theta depending on the precision
	  theta.mul(theta,eps);

	  
      	  FP_NR<FT> mu,mu_test;

      	  for (i=0; i<kappa; i++) {
	    
      	    mu.div(R.get(i,kappa),R.get(i,i));
      	    mu.abs(mu);

      	    mu_test.div(theta,R.get(i,i));
	    mu_test.abs(mu_test); // Jeu 12 mai 2016 12:55:13 CEST 
      	    mu_test.add(mu_test,one);
	    
      	    if (mu.cmp(mu_test) == 1) {
	      
      	      cout << " **** #tests = " << nblov << " **** Anomaly in size reduction, i = " << i << " ,kappa = " << kappa  << endl;
	      return -1;
      	    }
	    
      	  }
      	  somedone = 0;  // Here, should be size reduced, hence ok for continuing 
	  
      } // End test prec 
	    
      prev_nonstop = nonstop;
      
    }
      
    else 
      nonstop=0;

    
    // Heuristic test for not enough precision with respect to delta
    // Should be done much more efficiently
    
    // if ((nonstop==0) && (somedone ==1))  {
    
     
     
    // } // end test 

    
  } // end while 

  
  return somedone;

}


  /* -------------------------------------------------------------------------
   Size reduce sequential update of R and B   

   ------------------------------------------------------------------------- */

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline bool 
  Lattice<ZT,FT, MatrixZT, MatrixFT>::size_update(int kappa, int from_i, int to_i) { 

    bool somedone=false;
    
    int i;

    int nmax;
    
    Z_NR<ZT>  xz;
    
    long expo,lx;
  
    
    // Loop through the column 
    // -----------------------

    for (i=from_i; i>=to_i; i--){  
    
  
      x.div(R.get(i,kappa),R.get(i,i));
      x.rnd(x);
 
      
	
      if (x.sgn() !=0) {   // Non zero combination 
                           // --------------------
	lx = x.get_si_exp(expo);

		  
	nmax=structure[i]+1;
	
	// Cf fplll 
	// Long case 
	if (expo == 0) {

	  if (lx == 1) {

	    //compteur +=1;
	    somedone = 1;
	    
	    
	    R.subcol(kappa,i,i+1);
	    
	    B.subcol(kappa,i,nmax);
	    	    
	    if (transf) 
	      U.subcol(kappa,i,min(d,nmax));
	
	  } 
	  else if (lx == -1) {

	    //compteur +=1;
	    somedone = 1;
	   
	    
	    R.addcol(kappa,i,i+1);
	    
	    B.addcol(kappa,i,nmax);
		
	    if (transf) 
	      U.addcol(kappa,i,min(d,nmax));

	  } 
	  else { 

	    //compteur +=1;
	    somedone = 1;
	    
 
	    if (fast_long_flag == 1) {
	      
	      R.submulcol(kappa,i,x,i+1);
	      B.addmulcol_si(kappa,i,-lx,nmax);
	      if (transf)  
		U.addmulcol_si(kappa,i,-lx,min(d,nmax));
	      

	    } // end fast_long
	    else {
	      
	      set_f(xz,x);

	      R.submulcol(kappa,i,x,i+1);
	      
	      B.submulcol(kappa,i,xz,nmax);
	      if (transf)  
	     	U.submulcol(kappa,i,xz,min(d,nmax));
	    }

	    
	  } // end else expo ==0 and not 1 or -1
  
	} // end expo == 0 
	else {  // expo <> 0 

	  //compteur +=1;
	  somedone = 1;
	 
 
	  if (fast_long_flag == 1) {
	    
	    R.submulcol(kappa,i,x,i+1);
	    B.addmulcol_si_2exp(kappa,i,-lx,expo,nmax);
	    if (transf)  
	      U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));
	    
	    } // end fast_long
	  else {
	   
	    set_f(xz,x);  
	  
	    R.submulcol(kappa,i,x,i+1);	
	    B.submulcol(kappa,i,xz,nmax);
	    if (transf)  
	      U.submulcol(kappa,i,xz,min(d,nmax));
	  
	  } // end no long
	  
	} // end expo <> 0 

      } // Non zero combination 

    } // Loop through the column
    
    
    return somedone;

  }


/* -------------------------------------------------------------------------
   Size reduce update of R    

   ------------------------------------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline bool 
Lattice<ZT,FT, MatrixZT, MatrixFT>::size_update_R(vector<FP_NR<FT> >& vectx, int kappa, int from_i, int to_i) { 

  bool somedone=false;
    
  int i;
  
  long expo,lx;
  
  
  // Loop through the column 
  // -----------------------
  
  for (i=from_i; i>=to_i; i--){  
    
  
    x.div(R.get(i,kappa),R.get(i,i));
    x.rnd(x);
    
    vectx[i]=x;
    
    if (x.sgn() !=0) {     // Non zero combination 
                           // --------------------
      lx = x.get_si_exp(expo);
      
      
      // Cf fplll 
      // Long case 
      if (expo == 0) {
	
	if (lx == 1) {
	  
	  
	  somedone = 1;
	  
	  R.subcol(kappa,i,i+1);
	  
	} 
	else if (lx == -1) {
	  
	  
	  somedone = 1;
	  
	  R.addcol(kappa,i,i+1);
	  
	}
	else {

	  somedone = 1;
	  
	  R.submulcol(kappa,i,x,i+1); 
	  
	}
      }
      else { 

	somedone = 1;
	
	R.submulcol(kappa,i,x,i+1);
	
      } // end expo <> 0 
      
    } // Non zero combination 

  } // Loop through the column
    
    
  return somedone;

}


  
  
  
  
/* -------------------------------------------------------------------------
   Size reduce parallel update of B   

   ------------------------------------------------------------------------- */

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline bool 
  Lattice<ZT,FT, MatrixZT, MatrixFT>::psize_update_B(int kappa, int from_i, int to_i,  vector<FP_NR<FT> > vectx, int S) {

    
    int nmax;
  
    nmax=structure[from_i];
    
    for (int i=from_i+1; i>=to_i; i--)  
      nmax=max(nmax,structure[i]);

   
    nmax +=1;

#ifdef _OPENMP
#pragma omp parallel for 
#endif
    
    // Parallel loop
    // -------------
    
    for (int l=0; l<S; l++) { 

      int i;

      Z_NR<ZT>  xz;

      FP_NR<FT>  xf; // Old x 
    
      long expo,lx;
    
      int nloc;
      nloc=nmax/S;
      if ((l+1) <= nmax%S)
	nloc+=1;
      
      int ibeg;
      ibeg=(nmax/S)*l;
      if ((l+1) <= nmax%S)
	ibeg+=l;
      else
	ibeg+=nmax%S;

      
      // Loop through the column 
      // -----------------------

      for (i=from_i; i>=to_i; i--){  
    
  
	xf=vectx[i]; 

	
	if (xf.sgn() !=0) {   // Non zero combination 
	                      // --------------------
	  
	  lx = xf.get_si_exp(expo); // To put also through R 
		  
	  //nmax=structure[i]+1;
	
	  // Cf fplll 
	  // Long case 
	  if (expo == 0) {

	    if (lx == 1) {

	      B.psubcol(kappa,i,ibeg,nloc);
	    	    
	      // if (transf)   A FAIRE EN PARALLELE CF LONG 
	      // 	U.subcol(kappa,i,min(d,nmax));
	
	    } 
	    else if (lx == -1) {

	      B.paddcol(kappa,i,ibeg,nloc);
		
	      // if (transf) 
	      // 	U.addcol(kappa,i,min(d,nmax));

	    } 
	    else { 

	      if (fast_long_flag == 1) {
	      
		B.paddmulcol_si(kappa,i,-lx,ibeg,nloc);
		
		// if (transf)  
		//   U.addmulcol_si(kappa,i,-lx,min(d,nmax));
	      

	      } // end fast_long
	      else {
	      
		set_f(xz,xf);
	      
		B.psubmulcol(kappa,i,xz,ibeg,nloc);
		
		// if (transf)  
		//   U.submulcol(kappa,i,xz,min(d,nmax));
	      }

	    
	    } // end else expo ==0 and not 1 or -1
  
	  } // end expo == 0 
	  else {  // expo <> 0 


	    if (fast_long_flag == 1) {
	    
	      B.paddmulcol_si_2exp(kappa,i,-lx,expo,ibeg,nloc);
	      
	      // if (transf)  
	      // 	U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));
	    
	    } // end fast_long
	    else {
	   
	      set_f(xz,xf);  
	  
	      B.psubmulcol(kappa,i,xz,ibeg,nloc);
	      
	      // if (transf)  
	      // 	U.submulcol(kappa,i,xz,min(d,nmax));
	  
	    } // end no long
	  
	  } // end expo <> 0 

	} // Non zero combination 

      } // Loop through the column
    
    } // Parallel loop
    
    return true; // Somedone by R 
  }
  
  

  
/* ------------- */
/* Householder R */
/* ------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::householder_r(int kappa)
{

  nmaxkappa=structure[kappa]+1;

  int i,k,length; 

  // kappa == 0
  // ----------
  if (kappa==0) {

    if (col_kept[kappa]) {
     
      R.setcol(kappa,Bfp.getcol(kappa),nmaxkappa);
    
    }
     else {
      
       col_kept[kappa]=1; 

       Bfp.setcol(kappa,B.getcol(kappa),0,nmaxkappa);
       R.setcol(kappa,Bfp.getcol(kappa),nmaxkappa);
       fp_norm_sq(normB2[kappa], R.getcol(kappa), nmaxkappa);

     }
 
    kappamin[kappa]=kappa;

  }

  // kappa >0 
  // --------
  else {
    
    // Length of all the computations 
    
    length=nmaxkappa;
    
    // ---------------------------------------------
    // Re-use of previous orthononalization data 
    // ---------------------------------------------
    
    
    //if (descendu[kappa]>=1) {
      // Actually already implicitly considered in the next test 
      
    //}
    //else 
    if (col_kept[kappa]) { 
    
      // Keep the current norm and the floating point value of B 
      // Everything if not too big index decrease otherwise just a slice, for the index increase also 
      if  (((descendu[kappa] < 1) || (kappa-descendu[kappa] <= 0))) {

	// k=0
	if (kappamin[kappa]==0) {
	  
	  k=0;
	  scalarprod(VR(k,kappa), V.getcol(k,k), Bfp.getcol(kappa,k), length);
	  Rkept.fmasub(0,k,Bfp.getcol(kappa,k), V.getcol(k,k), VR(k,kappa), length); // k=0
	  
	  for (k=1; (k<kappa) && (k<nmaxkappa) ; k++) {
	    length--;
	    scalarprod(VR(k,kappa), V.getcol(k,k), Rkept.getcol(k-1,k), length);
	    Rkept.fmasub(k,k, Rkept.getcol(k-1,k), V.getcol(k,k), VR(k,kappa), length);  //  k-1 to k 
	  }
	}
	else {

	  k=0;
	  
	  Rkept.fmasub(k,k, Bfp.getcol(kappa,k), V.getcol(k,k), VR(k,kappa), length);  // k=0
	  
	  // (k < nmax kappa) added Mer 28 mai 2014 11:15:45 CEST for the rectangular case 
	  for (k=1; (k<kappamin[kappa]) && (k < nmaxkappa); k++)  { 
	    length--;
	    Rkept.fmasub(k,k, Rkept.getcol(k-1,k), V.getcol(k,k), VR(k,kappa), length);  // k-1 to k
	  } 
	  
	  for (k=kappamin[kappa]; (k<kappa)  && (k < nmaxkappa); k++) {
	    length--;
	    scalarprod(VR(k,kappa), V.getcol(k,k), Rkept.getcol(k-1,k), length);
	    Rkept.fmasub(k,k, Rkept.getcol(k-1,k), V.getcol(k,k), VR(k,kappa), length);  // k-1 to k 
	  }
	}
     
      } // endif not big decrease or increase 

    } // end colkept 
    // -----------------------------------------------------------
    // Complete re-computation  
    // -----------------------------------------------------------
    else { 

      
      col_kept[kappa]=1;  
      
      Bfp.setcol(kappa,B.getcol(kappa),0,nmaxkappa);
       
       
      fp_norm_sq(normB2[kappa], Bfp.getcol(kappa), nmaxkappa);
      
      // k =0 
      k=0;
      scalarprod(VR(k,kappa), V.getcol(k,k), Bfp.getcol(kappa,k), length);
      
      Rkept.fmasub(0,k, Bfp.getcol(kappa,k), V.getcol(k,k), VR(k,kappa), length); 
 
      // (k < nmax kappa) added Mer 28 mai 2014 11:15:45 CEST for the rectangular case
      for (k=1; (k<kappa) && (k < nmaxkappa); k++) {
	
	length--;
	
	scalarprod(VR(k,kappa), V.getcol(k,k), Rkept.getcol(k-1,k), length);
		
	Rkept.fmasub(k,k, Rkept.getcol(k-1,k), V.getcol(k,k), VR(k,kappa), length);  // de k-1 à k 
      }
      
      length = nmaxkappa; 

    } // endelse recomputation 
    // -----------------------------------------------------------
    
     
    // Dummy in the standard case, made special for the MatrixPE case 
   
    //for (i=0; i<kappa; i++) toR[i]=Rkept.get_non_normalized(i,i);
    // GV Mer 21 mai 2014 17:11:32 CEST for the rectangular case

    // Should be done more efficiently
    
    if (kappa < nmaxkappa) {
      
      for (i=0; i<kappa; i++) toR[i]=Rkept.get_non_normalized(i,i);
      for (i=kappa; i<nmaxkappa; i++) toR[i]=Rkept.get_non_normalized(i,kappa-1);
      
      R.setcol(kappa,&toR[0],nmaxkappa);
      
    }
    else {
      
      for (i=0; i< nmaxkappa; i++) toR[i]=Rkept.get_non_normalized(i,i);
    
      R.setcol(kappa,&toR[0],nmaxkappa);

      }
    
    kappamin[kappa]=kappa;

  } // else kappa !=0 

  
 return 0;
}

/* ------------- */
/* Householder V */
/* ------------- */
// nmaxkappa must be initialized e.g. through householder_r 

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::householder_v(int kappa) 
{

 int i;
  FP_NR<FT> s,norm,w,tmpdpe; 
  s=0;
  tmpdpe=0;
 
  //R.normalize(kappa,nmaxkappa);  // voir si nécessaire ? Rajouter en dummy si besoin aussi mpfr 

  w=R.get(kappa,kappa);


  if (w >=0) {

    fp_norm(s,R.getcol(kappa,kappa),nmaxkappa-kappa); 
    tmpdpe.neg(s);
    R.set(kappa,kappa,tmpdpe);  // On ne met pas à zéro, inutile, sauf pour getR
    
  }
  else {

    fp_norm(tmpdpe,R.getcol(kappa,kappa),nmaxkappa-kappa); // de la colonne 
    R.set(kappa,kappa,tmpdpe);
    
    s.neg(tmpdpe); 
  }

  w.add(w,s);

  s.mul(s,w);
  s.sqrt(s);

  V.div(kappa,kappa+1, R.getcol(kappa,kappa+1), s, nmaxkappa-kappa-1);
 
  // vraiment utile il n'y a pas deja des 0?   --> sans doute pas 12/04/11  
  for (i=nmaxkappa; i< n; i++) V.set(i,kappa,0.0);  // V à zéro car ré-utilisé plus loin ensuite (pas R); 


  tmpdpe.div(w,s);
    
  V.set(kappa,kappa,tmpdpe); 
  
  return 0; 
}


  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline unsigned int 
  Lattice<ZT,FT, MatrixZT, MatrixFT>::set_nblov_max(unsigned int nb) {

    nblov_max = nb;
    return nblov_max;

  } 

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline  int 
  Lattice<ZT,FT, MatrixZT, MatrixFT>::set_num_S(int S, int dthreshold) {

   
    num_S=S;
    par_reduce_threshold=dthreshold;
    
     
#ifdef _OPENMP
    omp_set_num_threads(num_S);
#endif
    
    return num_S;
  }
    
template<class ZT,class FT, class MatrixZT, class MatrixFT> inline unsigned int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::setprec(unsigned int prec) {

  // Re-initializations
  // Should be for each variable, non global? 
  // push and pop ? 
 
  unsigned oldprec;
  oldprec=getprec();

  mpfr_set_default_prec(prec);

  R.clear();
  R.resize(n,d);
  unsigned newprec;
  newprec=(R.get(0,0)).getprec(); 
  if (newprec == oldprec) cout << "Warning: in function setprec hlll, the change of precision has no effect" << endl; 

  Rkept.clear();
  Rkept.resize(n,d);

  V.clear();
  V.resize(n,d);

  VR.clear();
  VR.resize(d,d);

  Bfp.clear();
  Bfp.resize(n,d);

  normB2.clear();
  normB2.resize(d); 
  
  toR.clear();
  toR.resize(n); 
 
  // The old precision 
  return oldprec;
}


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline unsigned int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::getprec() {

  return (R.get(0,0)).getprec(); 
  
}


template<class ZT,class FT,class MatrixZT, class MatrixFT> inline ZZ_mat<ZT> Lattice<ZT,FT, MatrixZT, MatrixFT>::getbase()
{
   ZZ_mat<ZT> BB(n,d);
   for (int i=0; i<n; i++) 
     for (int j=0; j<d; j++) BB.Set(i,j,B(i,j)); // reprendre boucle sur les colonnes 

  
  return BB;
}


template<class ZT,class FT,class MatrixZT, class MatrixFT> inline MatrixZT Lattice<ZT,FT, MatrixZT, MatrixFT>::getmixedbase()
{
  return B;
}


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline ZZ_mat<ZT> Lattice<ZT,FT, MatrixZT, MatrixFT>::getU()
{
 
  if (transf) { 
    ZZ_mat<ZT> UU(d,d); 
    for (int i=0; i<d; i++) 
      for (int j=0; j<d; j++) UU.Set(i,j,U.get(i,j)); // reprendre boucle sur les colonnes 
    return UU;
  }
  else {
    cout << "*** Error, HLLL, the transformation matrix has not been computed" << endl;
    ZZ_mat<ZT> UU(0,0);
    return UU;
  }
}

  
template<class ZT,class FT, class MatrixZT, class MatrixFT> inline  matrix<FP_NR<FT> > Lattice<ZT,FT, MatrixZT, MatrixFT>::getR()
{
  matrix<FP_NR<FT> >  RR(d,d);
  FP_NR<FT> tmp;

  for (int i=0; i<min(n,d); i++) 
    for (int j=i; j<d; j++) {
      tmp=R.get(i,j);  // cf l'absence de const dans nr.cpp Set / Exp 
      RR.set(i,j,tmp); // reprendre boucle sur les colonnes 
      
    }
  for (int i=0; i<d; i++) 
    for (int j=0; j<i; j++) RR(i,j)=0.0;
  
  return RR;
}


// Constructeur 
// ------------


template<class ZT,class FT, class MatrixZT, class MatrixFT> void 
Lattice<ZT,FT, MatrixZT, MatrixFT>::init(int n, int d, bool forU) {

  int i,j;

  transf=forU;
  compteur=0;
  tmpcompt=0;

  tps_reduce=0;
  tps_householder=0;
  tps_prepare=0;
  tps_swap=0;
  nblov=0;
  nbswaps=0;
  tps_redB=0;

 
  R.resize(n,d);

  Rkept.resize(n,d);

  Bfp.resize(n,d);

  normB2.resize(d);
  toR.resize(n);
  
  V.resize(n,d);
  
  col_kept.resize(d+1); // +1 for the discovery test 
  descendu.resize(d);
  for (i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}

  kappamin.resize(d); // Lowest point down for kappa since last time
  for (j=0; j<d; j++) kappamin[j]=-1;

  VR.resize(d,d);

}


template<class ZT,class FT, class MatrixZT, class MatrixFT>
Lattice<ZT,FT, MatrixZT, MatrixFT>::Lattice(ZZ_mat<ZT> A, bool forU, int reduction_method, int long_flag) {

  
  n=A.getRows();
  d=A.getCols();

  init(n,d, forU); 

  int i,j;

  B.resize(n,d);  // Not in init for the mixed matrix case also 

  for (i=0; i<n; i++) 
    for (j=0; j<d; j++) 
      B(i,j)=A.Get(i,j);


  if (transf) {    // Not in init for the mixed matrix case also 
    
    U.resize(d,d);
    for (i=0; i<d; i++) U(i,i)=1;     
  }

  kappa_max = d;

  nblov_max = 4294967295;
  
  seysen_flag=reduction_method;

  fast_long_flag = long_flag;
    
  matrix_structure(structure, B, n,d);

  // For openmp, no threads by default 
  num_S=1;
  par_reduce_threshold=1000000;

 }


// With mixed matrix from mixed matrix 
// -----------------------------------

// template<class ZT,class FT, class MatrixZT, class MatrixFT>
// Lattice<ZT,FT, MatrixZT, MatrixFT>::Lattice(MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> > A, bool forU, int reduction_method) {

//     n=A.getRowsRT()+A.getRowsZT();
//     d=A.getCols();

//     init(n,d, forU); 
    
//     int i,j;

//     B.resize(A.getRowsRT(), A.getRowsZT(), d);

//     for (i=0; i<A.getRowsRT(); i++) 
//       for (j=0; j<d; j++) 
// 	B.setRT(i,j,A.getRT(i,j));

//     for (i=0; i<A.getRowsZT(); i++) 
//       for (j=0; j<d; j++) 
// 	B.setZT(i,j,A.getZT(i,j));
    
//     if (transf) {    // Not in init for the mixed matrix case also 
//       Z_NR<ZT> one;
//       one=1;
//       U.resize(0,d,d);
//       for (i=0; i<d; i++) U.set(i,i,one); 
//     }

//     nblov_max = 4294967295;
    

//     // CHECK this with base constructor changed 
//     seysen_flag=reduction_method;

//     if (reduction_method == DEF_REDUCTION) fast_long_flag = 1;
//     else if  (reduction_method == NO_LONG)  fast_long_flag = 0;
  
//     matrix_structure(structure, B, A.getRowsRT(), A.getRowsZT(), d);
// }


// ------------------------------------------------------------------
// Assigns a basis A 
// Assumes that the lattice has already been initialized 
// 
// -------------------------------------------------------------------  

  
template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
Lattice<ZT,FT, MatrixZT, MatrixFT>::assign(ZZ_mat<ZT> A) {

  nblov = 0;

  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) 
      B(i,j)=A.Get(i,j);
  
  matrix_structure(structure, B, n,d);
  
  for (int i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}
  for (int j=0; j<d; j++) kappamin[j]=-1;

  if (transf) {
    
    U.resize(d,d);
    for (int i=0; i<d; i++) U(i,i)=1; 
   
  }
  
}

template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
Lattice<ZT,FT, MatrixZT, MatrixFT>::assign(MatrixZT A) {

  nblov = 0;

  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) 
      B(i,j)=A.get(i,j);
  
  matrix_structure(structure, B, n,d);
  for (int i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}
  for (int j=0; j<d; j++) kappamin[j]=-1;

  if (transf) {
    
    U.resize(d,d);
    for (int i=0; i<d; i++) U(i,i)=1; 
   
  }
}


        
template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
Lattice<ZT,FT, MatrixZT, MatrixFT>::shift_assign(ZZ_mat<ZT> A, vector<int> shift, int sigma) {

  nblov = 0;

  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) 
      B(i,j).mul_2si(A(i,j),shift[i]);
  
  
  //    int ll = size_in_bits(B(0,0));
    
    //for (int i=0; i<n; i++)
    //for (int j=0; j<d; j++) 
    //	B(i,j).mul_2si(B(i,j),sigma - ll + 8);
    
  
  

  matrix_structure(structure, B, n,d);
  for (int i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}
  for (int j=0; j<d; j++) kappamin[j]=-1;

  if (transf) {
    
    U.resize(d,d);
    for (int i=0; i<d; i++) U(i,i)=1; 
   
  }
}


// ------------------------------------------------------------------
// Assigns a basis A with truncation 
// throw tau digits in the upper part: upperdim x d 
// keeps t digits globally 
//    if t=0 keeps evrything
//    if t<0 puts the identity in the lower part: (n-upperdim) x d 
//
// Assumes that the lattice has already been initialized 
// -------------------------------------------------------------------  

template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
Lattice<ZT,FT, MatrixZT, MatrixFT>::put(ZZ_mat<ZT> A, long upperdim, long t, long tau) {

  // Cf structure, colkept, kappamin 

  trunc<ZT, MatrixZT>(B, A, upperdim, d, n, t, tau);

  if (transf) {

    U.resize(d,d);
    for (int i=0; i<d; i++) U(i,i)=1; 

  }
}


template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
Lattice<ZT,FT, MatrixZT, MatrixFT>::mixed_put(MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> > A, long t, long tau) {

  // Cf structure, colkept, kappamin 

  mixed_trunc<matrix, FP_NR<mpfr_t>, Z_NR<ZT> >(B, A, t, tau);

  if (transf) {
    
    Z_NR<ZT> one;
    one=1;
    U.resize(0,d,d);
    for (int i=0; i<d; i++) U.set(i,i,one); 
  }
  
}


// After initialization 
// Multiply the first m rows by 2^sigma 

template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
Lattice<ZT,FT, MatrixZT, MatrixFT>::shift(ZZ_mat<ZT> A, long m, long lsigma) {
  // Cf structure, colkept, kappamin 
  int i,j;

  for (i=0; i<m; i++)
    for (j=0; j<d; j++)
      B(i,j).mul_2exp(A(i,j), lsigma); 

  if (transf) {

    U.resize(d,d);
    for (int i=0; i<d; i++) U(i,i)=1; 

  }
}

template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
Lattice<ZT,FT, MatrixZT, MatrixFT>::shiftRT(long lsigma) {
  // Cf structure, colkept, kappamin 
  B.shift(lsigma); 

  if (transf) {

    Z_NR<ZT> one;
    one=1;
    U.resize(0,d,d);
    for (int i=0; i<d; i++) U.set(i,i,one); 

  }
}



//***************************************************************************
// Reduction test 
// 
// Relevant iff mpfr is used (for the precision change)  
//  
// Matrix interface () deprecated, will note work with -exp 
//

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void Lattice<ZT,FT, MatrixZT, MatrixFT>::isreduced(double deltain) {

  cout << endl << "Testing reduction" << endl;
  
  // We launch HLLL, "reduced" if no swap
  // ************************************

  unsigned int oldprec,newprec;
  int k;
  FP_NR<FT> s1,s2,s,delta; 
  Z_NR<ZT> nb1;

  oldprec=getprec();

  if (d<=20) setprec(53);  
  else setprec(2*d);  // ******* TO TUNE  

  newprec=getprec();

  if (newprec == oldprec) cout << "Warning: in function isreduced, the change of precision has no effect" << endl; 

  // Premier calcul avant test et further éventuel 
  // ---------------------------------------------

  householder();

  s1.mul(R(0,1),R(0,1));
  s2.mul(R(1,1),R(1,1));
  s.mul(R(0,0),R(0,0));

  s1.add(s1,s2);
  delta.div(s1,s);
  if (delta >=1) delta=1.0; 

  for (k=1; k<d-1; k++) {

    s1.mul(R(k,k+1),R(k,k+1));
    s2.mul(R(k+1,k+1),R(k+1,k+1));
    s.mul(R(k,k),R(k,k));
    s1.add(s1,s2);
    s.div(s1,s);

    if (s < delta) delta=s;
  }
  
  cout << "n = "<< n << endl;
  cout << "d = "<< d << endl;
  cout << "(Householder) delta is about "; 
  delta.print(); 
  cout << endl;

  fp_norm_sq(nb1,B.getcol(0),n);
  set_z(s,nb1);
  s.sqrt(s);
      
  cout << "(Householder) ||b_1|| is ";
  s.print();
  cout  << endl; 

  s=R(0,0);
  for (k=1; k<d; k++) 
    s.mul(s,R(k,k));

  s.abs(s);
  cout << "(Householder) Vol(L) is about ";
  s.print();
  cout  << endl; 

  FP_NR<FT> eta,theta,alpha;  // Pas mpfr 

  theta=0.00001;

  eta=0.0;

  int i,j;

  FP_NR<FT> v,w; // pas mpfr 

  for (i=0; i<d-1; i++)
    for (j=i+1; j<d; j++) {
      v.abs(R(i,j));
      w.mul(theta,R(j,j));
      v.sub(v,w);
      if (v >0) {
	v.div(v,R(i,i));
	if (v > eta) eta=v;
      }
    }

  alpha.mul(theta,theta);
  v=1.0;
  alpha.add(alpha,v);

  alpha.mul(alpha,delta);

  v.mul(eta,eta);
  alpha.sub(alpha,v);
  alpha.sqrt(alpha);

  v.mul(eta,theta);
  alpha.add(v,alpha);

  v.mul(eta,eta);
  v.sub(delta,v);

  alpha.div(alpha,v);


  cout << "eta : ";
  eta.print();
  cout << "     theta : ";
  theta.print();
  cout << "    alpha : "; 
  alpha.print(); 
  cout << endl; 

  // Reduction  
  // ---------

  nblov=0;

  hlll(deltain);

  if (nblov== ((unsigned int) d-1)) cout << "Seems to be reduced (up to size-reduction)" << endl << endl;

  if (nblov != ((unsigned int) d-1)) {
    cout << endl;
    cout << "!! Does not seem to be reduced" << endl;
    cout << "     #non-trivial swaps is " << nblov-d+1 << " > " << 0 << endl << endl; 
  

    s1.mul(R(0,1),R(0,1));
    s2.mul(R(1,1),R(1,1));
    s.mul(R(0,0),R(0,0));

    s1.add(s1,s2);
    delta.div(s1,s);
    if (delta >=1) delta=1.0; 

    for (k=1; k<d-1; k++) {

      s1.mul(R(k,k+1),R(k,k+1));
      s2.mul(R(k+1,k+1),R(k+1,k+1));
      s.mul(R(k,k),R(k,k));

      s1.add(s1,s2);
      s.div(s1,s);
    
      if (s < delta) delta=s;
  
    }

    cout << "(further with " << nblov-d+1 << " tests) " << "delta is about "; 
    delta.print(); 
    cout << endl;


    fp_norm_sq(nb1,B.getcol(0),n);
    set_z(s,nb1);
    s.sqrt(s);
      
    cout << "(further with " << nblov-d+1 << " tests) " << "||b_1|| is ";
    s.print();
    cout  << endl; 


    s=R(0,0);
    for (k=1; k<d; k++) 
      s.mul(s,R(k,k));

    s.abs(s);
    cout << "(further with " << nblov-d+1 << " tests) " << "Vol(L) is about ";
    s.print();
    cout  << endl; 

    

    theta=0.001;
    
    eta=0.0;
  
  for (i=0; i<d-1; i++)
    for (j=i+1; j<d; j++) {
      v.abs(R(i,j));
      w.mul(theta,R(j,j));
      v.sub(v,w);
      v.div(v,R(i,i));
      if (v > eta) eta=v;
    }

  alpha.mul(theta,theta);
  v=1.0;
  alpha.add(alpha,v);

  alpha.mul(alpha,delta);

  v.mul(eta,eta);
  alpha.sub(alpha,v);
  alpha.sqrt(alpha);

  v.mul(eta,theta);
  alpha.add(v,alpha);

  v.mul(eta,eta);
  v.sub(delta,v);

  alpha.div(alpha,v);


  cout << "eta : " << eta << "    theta : " << theta << "    alpha : " << alpha << endl << endl; 

  } // Was not reduced 

  setprec(oldprec);
}


//***************************************************************************
// Upper bound on log[2] Condition number || |R| |R^-1| || _ F
// 
// R is not assumed to be known 
//
//*************************************************************************** 
// Matrix interface () deprecated, will note work with -exp 
//

//#define DEFAULT_PREC 0 
//#define UNKNOWN_PREC 1
//#define CHECK 1
// PREC IF >=2 

//#define TRIANGULAR_PROPER 1   
//#define ANY 0 

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline long  
Lattice<ZT,FT, MatrixZT, MatrixFT>::lcond(int tproper, int flagprec,  int flagheur) {


  if (flagprec == DEFAULT_PREC) {

    // Default prec and triangular 
    // ***************************

    if (tproper == TRIANGULAR_PROPER) {

      int i,j;
      
      for (j=0; j<d; j++)
	R.setcol(j,B.getcol(j),0,j+1); 

      FP_NR<FT>  tmp,ttmp;

      // Size-reducedness 
  
      FP_NR<FT>  eta;
      eta=0.0;

      for  (j=1; j<d; j++) 

	for (i=0; i<j; i++) {
	  tmp.div(R(i,j),R(i,i));
	  tmp.abs(tmp);
	  if (tmp.cmp(eta) > 0) eta=tmp;
	}

      // diag quo
        
      FP_NR<FT>  maxquo;
      maxquo=0.0;

      for  (j=0; j<d; j++) 
	for (i=0; i<d; i++) {
	  tmp.div(R(j,j),R(i,i));
	  tmp.abs(tmp);
	  if (tmp.cmp(maxquo) > 0) maxquo=tmp;
	}
  
      // Approx cond all together 
     
      long condb; 
 
      condb=maxquo.exponent();
      tmp=1.0;
      eta.add(tmp,eta);
      condb += (d-1) + eta.exponent();

      Z_NR<long> dd;
      dd=d;
      condb += 1 + dd.exponent();

      return condb;


    } // end if triangular and viewed as proper 

    // Defaut prec with no structure 
    // *****************************

    else { // No structure but default prec  

      int i,j,k;

      // Direct computation for R with householder 
      householder();

      // Do not change for the template class MatrixZT, class MatrixFT that is different e.g. double and mpfr 
      MatrixFT aR;
      aR.resize(d,d);
      
      for (i=0; i<d; i++)
	for (j=0; j<d; j++) 
	  aR(i,j).abs(R(i,j));

      // Inversion 
      // ---------

      MatrixFT iR;
      iR.resize(d,d);

      for (i=0; i<d; i++) iR(i,i)=1.0; 

      // Higham p263 Method 1

      FP_NR<FT> t,one;
      one=1.0;

      for (j=0; j<d; j++) {
	  iR(j,j).div(one,aR(j,j));
	  for (k=j+1; k<d; k++) {
	    t.mul(iR(j,j),R(j,k));
	    iR(j,k).neg(t);
	  }
	  
	  for (k=j+1; k<d; k++) {
	    for (i=j+1; i<k; i++) {
	      t.mul(iR(j,i),R(i,k)); 
	      iR(j,k).sub(iR(j,k),t); 
            }
	    iR(j,k).div(iR(j,k),R(k,k)); 
	  }
      }

      // Abs iR --> iR 
      for (i=0; i<d; i++)
	for (j=0; j<d; j++) 
	  iR(i,j).abs(iR(i,j));

      // aR * iR 
      
      MatrixFT prod;
      prod.resize(d,d);
      
      for (i=0; i<d; i++)
	for (j=0; j<d; j++) {
	  prod(i,j).mul(aR(i,0),iR(0,j));
	  for (k=1; k<d ; k++) 
	    prod(i,j).addmul(aR(i,k),iR(k,j));
	}

      // Norm 

      FP_NR<FT>  cc,c2;
      cc=0.0; 
      
      for (i=0; i<d; i++)
	for (j=0; j<d; j++) 
	  cc.addmul(prod(i,j),prod(i,j));

      cc.sqrt(cc); 

      return(cc.exponent());
    } 
    
  } // end if default_prec 
  
  // Using a user defined prec, no structure 
  // ***************************************
  // Only through mpfr here 

  else if ((flagprec >=2) && (flagheur ==0)) {
    
    unsigned oldprec;
    oldprec=getprec();

    mpfr_set_default_prec(flagprec);

    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(getbase());

    return(L.lcond(ANY, DEFAULT_PREC));

    mpfr_set_default_prec(oldprec);

  } // and else using the user prec, no check 
  
  // Heuristic check by increasing the precision  
  // *******************************************
  // To tune the speed of precision increase 

  else if ((flagprec >= 2) && (flagheur == CHECK)) {

    unsigned oldprec,prec=0;

    oldprec=getprec();

    bool found = false;

    int cond,oldcond=0;

   

    while (found == false) {
      
      prec += flagprec;

 cout << "prec = " << prec << endl;
      mpfr_set_default_prec(prec);

      Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(getbase());
      
      cond = L.lcond(ANY, DEFAULT_PREC);

      if (cond == oldcond) 
	found = true;
      else 
	oldcond = cond;
    }
    
    mpfr_set_default_prec(oldprec);

    return cond; 

  }  // and else using the user prec, with check

  // Unknown prec 
  // ************
  // To tune  
  else if (flagprec == UNKNOWN_PREC) {

    long bits = maxbitsize(getbase());

    unsigned oldprec;
    oldprec=getprec();
    
    mpfr_set_default_prec(2*d*bits);
    
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(getbase());
    
    return(L.lcond(ANY, DEFAULT_PREC));
    
    mpfr_set_default_prec(oldprec);


  } // and unknown 

  return 0; //cc
}




/* --------------------------------------------- */
/* Householder complet */
/* --------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::householder()
{

  int i,k,kappa;
  FP_NR<FT> nrtmp,s,w; 
  
  nrtmp=0;
  s=0;

  

    for (kappa=0; kappa<d; kappa++) {

      R.setcol(kappa,B.getcol(kappa),0,n);

      for (k=0; k<kappa; k++) {
	scalarprod(nrtmp, V.getcol(k,k), R.getcol(kappa,k), n-k);
	R.fmasub(kappa,k,R.getcol(kappa,k), V.getcol(k,k), nrtmp, n-k); 
      }


      w=R.get(kappa,kappa);

      if (w >=0) {
	fp_norm(s,R.getcol(kappa,kappa),n-kappa); 
	nrtmp.neg(s);
	R.set(kappa,kappa,nrtmp);    
      }
      else {
	fp_norm(nrtmp,R.getcol(kappa,kappa),n-kappa); // de la colonne
	R.set(kappa,kappa,nrtmp);
	s.neg(nrtmp);  
      }

      w.add(w,s);
      s.mul(s,w);
      s.sqrt(s);

      V.div(kappa,kappa+1, R.getcol(kappa,kappa+1), s, n-kappa-1);

      nrtmp.div(w,s);
      V.set(kappa,kappa,nrtmp); 

      for(i=kappa+1; i<d; i++)  R.set(i,kappa,0.0); 
       
    }  // sur kappa 
   
    return 0; 
}


  
  template<class ZT,class FT, class MatrixZT, class MatrixFT>  int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::hlllp(double delta, int dthreshold, int compt) { 

  
  int kappa=1,i;
  int prevkappa=-1; // For the looping test between tow indices 
  vector<FP_NR<FT> >  prevR(d);

  Timer cpu;
  cpu.clear();

  Timer tsize;
  tsize.clear();

  FP_NR<FT> newt; //testaccu;

  FP_NR<FT> deltab,lovtest;
  deltab=delta;   // TO SEE 

  FP_NR<FT> tmpswap;

  FP_NR<FT> s,sn; // newt test 

  int flag_reduce=0; // No convergence in reduce 


  for (i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}

  // **********
  // int kmax=2;

  // Timer time;

  // time.start();

  // int nb=0;
  
  //************
  
  cpu.start();

  while ((kappa < d) && (nblov < nblov_max)) 
    {

      	  
      if (((nblov%800000)==0) && (nblov > 0))   cout << nblov << " tests" << endl; 
      
      if (kappa == 1) { 
	for (i=0; i<d; i++) kappamin[i]=min(kappamin[i],0);
	householder_r(0);  
	householder_v(0);
	
      } 
      else for (i=kappa+1; i<d; i++) kappamin[i]=min(kappamin[i],kappa); // lowest visit after kappa 

      Timer debg;       
      debg.clear();

     
      
      debg.start();
      if (seysen_flag < 1) 
      	flag_reduce=hsizereduce(kappa);
      else {
	if (kappa > dthreshold) flag_reduce=newseysenreduce(kappa,dthreshold,compt);
	else  flag_reduce=seysenreduce(kappa);
      }
      debg.stop();
      if (kappa > dthreshold) tsize +=debg;
     
            
      if (flag_reduce==-1) return(-1);
      

      //newt=normB2[kappa];  // The newt test was like that in old non exp 
      //for (i=0; i<=kappa-2; i++) newt.submul(R.get(i,kappa),R.get(i,kappa));

      lovtest.mul(R.get(kappa-1,kappa-1),R.get(kappa-1,kappa-1));
      lovtest.mul(deltab,lovtest);

      nblov+=1;
    
      fp_norm(s,R.getcol(kappa,kappa),structure[kappa]+1-kappa);
      
      s.mul(s,s);

      
      
      sn.mul(R.get(kappa-1,kappa),R.get(kappa-1,kappa));
      newt.add(s,sn);
 

      // ****************
      //   UP    UP    UP 
      // ****************
    
      if (lovtest <= newt) {

	
	// ICI TESTS
	// if (((kappa) < d) && (col_kept[kappa+1] == 0)) {
	//   cpu.stop();
	//   cout << "Discovering vector " << kappa +1 << "/" << d << "   " << cpu << endl;
	//   cpu.start(); 
	// } 
							      
	
	householder_v(kappa);   // The first part of the orthogonalization is available 

	
	// Heuristique precision check : when R(kappa-1,kappa-1) increases in a 2x2 up and down  
	// ------------------------------------------------------------------------------------
	if (prevkappa==kappa+1) {  
	  FP_NR<FT> t;
	  t.abs(R.get(kappa,kappa));
       
	  if (t > prevR[kappa]) {
	
	    cout << " **** #tests = " << nblov << " **** Anomaly: the norm increases for kappa = " << kappa << endl;
	 
	    return -1;
	  }

	}

	

	descendu[kappa]=0;
	if (kappa==1) descendu[0]=0;
	
	prevkappa=kappa; 
	prevR[kappa].abs(R.get(kappa,kappa)); 
	
	kappa+=1; 



      } // End up 

      // ****************
      //   DOWN   DOWN 
      // ****************

      else {

	if (kappa==1) descendu[0]=0 ;
	else descendu[kappa-1]=descendu[kappa]+1;
	descendu[kappa]=0;

	nbswaps+=1;
       
	B.colswap(kappa-1,kappa);
	
	if (transf) U.colswap(kappa-1,kappa);
	
	Bfp.colswap(kappa-1,kappa);

	structure[kappa-1]=structure[kappa];

	VR.colswap(kappa-1,kappa);
 
	tmpswap=normB2[kappa];  normB2[kappa]=normB2[kappa-1]; normB2[kappa-1]=tmpswap;

	prevkappa=kappa; 
	kappa=max(kappa-1,1);  

      }

    } // End main LLL loop 

  // DBG
  cout << "--------------------------------------------------" << endl; 
  cout << endl << "*** tsize: " << endl;
  tsize.print(cout);
  cout << endl;
  cout << "--------------------------------------------------" << endl; 
   
  return 0;
  
};


/* -------------------------------------------------------------------------
   Seysen size reduction 

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm

   EXPERIMENTAL SEYSEN 

   ------------------------------------------------------------------------- */


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::newseysenreduce(int kappa, int dthreshold, int compt) { 

  //OMPTimer dbg;
  //dbg.clear();
  
  FP_NR<FT>  qq;
  
  vector<FP_NR<FT> > vectx(kappa);  
  
  vector<FP_NR<FT> > tmpcolR(kappa);  
 
  MatrixFT tR; 
  tR.resize(kappa,2); 

  int i,w=0;

  bool nonstop=1;
  bool somedone=0;

  vector<bool> bounded(kappa);
  
  int restdim=0; // Remaining dimension after the current block 

  
  householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 

  int bdim,ld,tdig,indexdec;
  
  while (nonstop) {  // LOOP COLUMN CONVERGENCE

    
    w++;

    somedone = 0;

    ld=1; indexdec=0; // Décalage d'indice
 

    while (ld <=kappa) {
      
      tdig=(kappa/ld)%2;
      
      if (tdig==0) bdim =0; else bdim = ld;

      // -----------------------------------------------------
      // Boucle sur la partie de la colonne correspond au bloc 
      // -----------------------------------------------------

      // vectxz rounding of vectx
      //  column -  (prev col) * vectxz (xz ré-utilisé) 
      // On peut travailler sur place en remontant dans la colonne kappa de R 

      // On calcule vectx et on arrondit au fur et à mesure

      restdim=kappa-indexdec-bdim;

      for (i=kappa-1-indexdec; i>=restdim; i--) 
	tmpcolR[i]=R.get(i,kappa);

      // À mettre Bien meilleur en dpe 
      //tR.setcol(1,R.getcol(kappa,restdim), bdim);

      compteur =0;
      
      for (i=kappa-1-indexdec; i>=restdim; i--){
	
	vectx[i].div(tmpcolR[i],R.get(i,i));

	// À mettre dpe 
	//vectx[i].div(tR.get(i-restdim,1),R.get(i,i));
	
	qq.abs(vectx[i]);
	if (qq.cmp(0.501) == 1) {
	  bounded[i]=0;
	  somedone = 1;
	  compteur +=1;
	}
	else bounded[i]=1;

	// Faire une opération vectorielle 
	for (int k=restdim; k<i; k++) 
	  tmpcolR[k].submul(R.get(k,i),vectx[i]);

	// À mettre dpe
	//tR.setcol(0,R.getcol(i,restdim), i-restdim); // i-restdim longueur ? 
	//tR.submulcol(1,0,vectx[i],i-restdim);	
	 
      } // end calcul de la transfo 
       

      // Et on applique la transformation  
      // --------------------------------

      int nthreads=1;

     
      if ((restdim > dthreshold) && (compteur > compt))  
	nthreads=2;
      
      //dbg.start();
      //#pragma omp parallel sections num_threads(nthreads)
      {

	//#pragma omp section
     
	{ // Update of R
	  // -----------

	  FP_NR<FT> x;
	  Z_NR<ZT>  xz;
	  long expo,lx;
	
	  for (int i=kappa-1-indexdec; i>= restdim; i--){

	    x.rnd(vectx[i]);
	 
	    if (bounded[i]==0) {
	  
	      lx = x.get_si_exp(expo);
	      
	      // Cf fplll 
	      // Long case 
	      if (expo == 0) {
	    
		if (lx == 1) {

		  R.subcol(kappa,i,restdim);
	      
		} 
		else if (lx == -1) {

		  R.addcol(kappa,i,restdim);
		  
		} 
		else { 

		  if (fast_long_flag == 1) {
	      
		    R.submulcol(kappa,i,x,restdim);
		  		
		  } // end fast_long
		  else {
	      
		    R.submulcol(kappa,i,x,restdim);	
		  }	      
		} 
	    
	      } // end expo == 0 
	      else {  // expo <> 0 

	  
		if (fast_long_flag == 1) {

		  R.submulcol(kappa,i,x,restdim);
		 
	    
		} // end fast_long
		else {
	   
		  set_f(xz,x);  
	  
		  R.submulcol(kappa,i,x,restdim);	
		  
	  
		} // end no long
	    
	      } // end expo <> 0
	      
	    } // Non zero combination, if on the column
	    
	  } // end application de la transformation, i loop on all the columns 

	} // End pragma update of R 
	
	//#pragma omp section
     
	{ // Update of B
	  // -----------

	  FP_NR<FT> x;
	  Z_NR<ZT>  xz;
	  long expo,lx;
	  int nmax; // De la structure triangulaire
	
	  for (int i=kappa-1-indexdec; i>= restdim; i--){

	    x.rnd(vectx[i]);

	 
	    if (bounded[i]==0) {
	  
	      lx = x.get_si_exp(expo);
	      
	      nmax=structure[i]+1;
	  
	      // Cf fplll 
	      // Long case 
	      if (expo == 0) {
	    
		if (lx == 1) {

		  // R.subcol(kappa,i,restdim);
	      
		  B.subcol(kappa,i,nmax);
	      
		  if (transf) 
		    U.subcol(kappa,i,min(d,nmax));
	      
		} 
		else if (lx == -1) {

		  // R.addcol(kappa,i,restdim);
	      
		  B.addcol(kappa,i,nmax);
	      
		  if (transf) 
		    U.addcol(kappa,i,min(d,nmax));

		} 
		else { 

		  if (fast_long_flag == 1) {
	      
		    // R.submulcol(kappa,i,x,restdim);
		    B.addmulcol_si(kappa,i,-lx,nmax);
		    if (transf)  
		      U.addmulcol_si(kappa,i,-lx,min(d,nmax));
		
		  } // end fast_long
		  else {
		
		    set_f(xz,x);  
	      
		    //  R.submulcol(kappa,i,x,restdim);	
		    B.submulcol(kappa,i,xz,nmax);
		    if (transf)  
		      U.submulcol(kappa,i,xz,min(d,nmax));
		  }	      
		} 
	    
	      } // end expo == 0 
	      else {  // expo <> 0 

	  
		if (fast_long_flag == 1) {

		  // R.submulcol(kappa,i,x,restdim);
		  B.addmulcol_si_2exp(kappa,i,-lx,expo,nmax);
		  if (transf)  
		    U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));
	    
		} // end fast_long
		else {
	   
		  set_f(xz,x);  
	  
		  // R.submulcol(kappa,i,x,restdim);	
		  B.submulcol(kappa,i,xz,nmax);
		  if (transf)  
		    U.submulcol(kappa,i,xz,min(d,nmax));
	  
		} // end no long
	    
	      } // end expo <> 0
	      
	    } // Non zero combination, if on the column
	    
	  } // end application de la transformation, i loop on all the columns 


	  
	} // End pragma update of B 

     
      } // End pragma parallel application of transformation

      //dbg.stop();
      if (nthreads ==2) {
	cout << "-----------------" << endl;
	//cout << dbg << endl;
	
      }

      indexdec+=bdim;     
      ld=ld*2;

    } // End loop on log blocks 

    
    if (somedone) {
       
            
      col_kept[kappa]=0;

      householder_r(kappa);  
     
    }
    else 
      nonstop=0;
     

  } // end while 
  
  return somedone;

};




} // end namespace hplll

#endif 

