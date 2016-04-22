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


#ifndef HPLLL_SLLL_CC
#define HPLLL_SLLL_CC

#ifdef _OPENMP
#include <omp.h>
#endif 

// *******************************************************************************
// 
//  S >= 2 parallel segments, hence nb_blocks = 2*S blocks
//
//  Divisibility of the dimension by nb_blocks in ensured by the constructor
//     hence S as the parameter of hlll should divide S given to the constructor
//
//  nbthreads should divide S (??) 
//
// *******************************************************************************


namespace hplll {   

template<class ZT,class FT, class MatrixZT, class MatrixFT>  int 
SLattice<ZT,FT, MatrixZT, MatrixFT>::hlll(double delta, int S, int nbthreads, unsigned int lovmax) { 

  int i,j,k;

#ifdef _OPENMP
  OMPTimer time,totime;
  OMPTimer redtime,eventime,oddtime,prodtime1,prodtime2,esizetime,osizetime;
  
  omp_set_num_threads(nbthreads);
#else 
  Timer time,totime;
  Timer redtime,eventime,oddtime,prodtime1,prodtime2,esizetime,osizetime;
#endif

  time.clear();
  redtime.clear();
  eventime.clear();
  oddtime.clear();
  prodtime1.clear();
  prodtime2.clear();
  esizetime.clear();
  osizetime.clear();
  totime.clear();
 
    
  totime.start();
  
  int nb_blocks=2*S;
  
  int bdim;    // Dimension of each block 
  // Assume that d is a multiple of nb_blocks  
  
  bdim = d/nb_blocks;
  
  int iter;
  
  bool stop=0;

  int condbits=53;
     
  // ************************************
  //   
  //   Main odd-even phases loop on iter
  //
  // ************************************

  // Computed through size reduction afterwards
  householder(); 

  for (iter=0; stop==0; iter++) {
  //for (iter=0; iter < 1; iter++) {

	    
    stop=1;
       

    // col_kept ??

    
    // The integer block lattice: truncation of the floating point R
    // 0 triangulaire ici car householder global

 
    set_f(RZ,R,condbits);
 
    
    for (i=1; i<d; i++)
      for (j=0;j<i;j++)
	RZ(i,j)=0;

    
    
    // In case one column is badly conditioned 
    // Make this clean, how?
 
    for (i=0; i<d; i++)
      if (RZ(i,i).sgn() ==0) RZ(i,i)=1;

   
    
    // Even reductions
    // ***************

    setId(U_even); // Pas nécessaire pour even sans doute 
    
    //cout << endl <<  "--- Even reductions" << endl;
    time.start();
    
#ifdef _OPENMP
#pragma omp parallel for 
#endif

    for (k=0; k<S; k++) { 
      {
	
	Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  BR(getblock(RZ,k,k,S,0),TRANSFORM,DEF_REDUCTION);
	
	BR.set_nblov_max(lovmax);
	BR.hlll(delta);
	
	//cout << endl << "Swaps (" << 2*k << "): " << BR.nbswaps << endl; 
	swapstab[2*k]+=BR.nbswaps;
	nbswaps+=BR.nbswaps;
	
	putblock(U_even,BR.getU(),k,k,S,0);

	//cout << "Nb threads: " << omp_get_num_threads() << endl; 

      }
    }

#ifdef _OPENMP
#pragma omp barrier 
#endif
    
    time.stop();
    redtime+=time;
    eventime+=time; 

    
    // Update after the block even reductions
    // **************************************

    stop= (stop &&  isId(U_even));

        
    time.start();

    
    pmatprod(S,0); 

     
    time.stop();
    prodtime1+=time;
    
    // Even size reduction
    // *******************
    
    time.start();
    
    for (i=0; i<d; i++) {
      col_kept[i]=0;
      descendu[i]=0;
    }
    
    householder_r(0);
    householder_v(0);
    
    for (i=1; i<d; i++) {
      
      hsizereduce(i);
      householder_v(i);

    }

    time.stop();
    esizetime+=time;

    
    // Odd reductions
    // **************

    setId(U_odd);
    
    set_f(RZ,R,condbits);

    // Forcer les zéros si pas mis dans householder_v par ex 
    for (i=1; i<d; i++)
      for (j=0;j<i;j++)
	RZ(i,j)=0;
    
    // DBG
    for (i=0; i<d; i++)
      if (RZ(i,i).sgn() ==0) RZ(i,i)=1;
    
    time.start();

    //cout << endl <<  "--- Odd reductions" << endl;
    
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    
    for (k=0; k<S-1; k++) {

      
      Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > BR(getblock(RZ,k,k,S,bdim),TRANSFORM,DEF_REDUCTION);
      
      BR.set_nblov_max(lovmax);
      BR.hlll(delta);
      //cout << endl << "Swaps (" << 2*k+1 << "): " << BR.nbswaps << endl;
      swapstab[2*k+1]+=BR.nbswaps;
      nbswaps+=BR.nbswaps;
      
      putblock(U_odd,BR.getU(),k,k,S,bdim);	   
      	  
    }
       
    
    time.stop();
    redtime+=time; 
    oddtime+=time; 
    
    // Update after the block odd reductions
    // *************************************

    stop = (stop &&  isId(U_odd));

    time.start();
    
    //matprod_in(B,U_odd);
    pmatprod(S,1);
    
    time.stop();
    prodtime2+=time;
    
    
    // Odd size reduction
    // ******************
    
    time.start();
    
    for (i=0; i<d; i++) {
      col_kept[i]=0;
      descendu[i]=0;
    }
    
    householder_r(0);
    householder_v(0);
    
    for (i=1; i<d; i++) {
      
      hsizereduce(i);
      householder_v(i);

    }

    time.stop();
    osizetime+=time;
    
  } // End odd-even iter until reduced  


  totime.stop();
  
  // cout << endl;
  // //cout << " Householder: " << qrtime << endl;
  // //cout << " Re-ortho: " << orthotime  << endl;
  // cout << " Reductions: " << redtime << endl;
  // cout << "   Even reductions: " << eventime << endl;
  // cout << "   Odd reductions: " << oddtime << endl;
  // cout << " Products: " << prodtime1 << endl;
  // cout << "           " << prodtime2 << endl;
  // //cout << "           " << prodtime3 << endl;
  // //cout << "           " << prodtime4 << endl;
  // cout << " Even size reds: " << esizetime  << endl;
  // cout << " Odd size reds: " << osizetime  << endl;
  // cout << " Total time:  " << totime << endl;  
  // //cout << endl << " Special chrono:" << special << endl << endl;
  // cout << endl << "Swaps: " << swapstab << endl;
  
      
  return 0;
  
};


/* -------------------------------------------------------------------------
   Seysen size reduction 

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm

   ------------------------------------------------------------------------- */


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
SLattice<ZT,FT, MatrixZT, MatrixFT>::seysenreduce(int kappa) { 

      
  nmaxkappa=structure[kappa]+1;

  FP_NR<FT> approx;
  approx=0.01;
   
  FP_NR<FT>  qq;
  
  FP_NR<FT> x,t,tmpfp;
  Z_NR<ZT>  xz,tmpz;

  long expo,lx;

  vector<FP_NR<FT> > vectx(kappa);  
  
  vector<FP_NR<FT> > tmpcolR(kappa);  
 
  int i,k,w=0;

  bool nonstop=1;
  bool somedone=0;

  vector<bool> bounded(kappa);
  
  int restdim=0; // Remaining dimension after the current block 

  int nmax; // De la structure triangulaire 


  //int whilemax=10000;  // Convergence problem if too low precision  

  // To see / prec problem 
  //col_kept[kappa]=0;
  
  householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 

  int bdim,ld,tdig,indexdec;
  
  while (nonstop) {  // LOOP COLUMN CONVERGENCE

    
    w++;

    somedone = 0;

    //compteur += 1;
    
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

	
      for (i=kappa-1-indexdec; i>=restdim; i--){
	
	vectx[i].div(tmpcolR[i],R.get(i,i));
	
 
	qq.abs(vectx[i]);
	if (qq.cmp(0.501) == 1) bounded[i]=0; else bounded[i]=1;

	
	
	// Faire une opération vectorielle 
	for (k=restdim; k<i; k++) 
	  tmpcolR[k].submul(R.get(k,i),vectx[i]);
	
	 
      } // end calcul de la transfo 
       

      // Et on applique la transformation  
      // --------------------------------
      for (i=kappa-1-indexdec; i>= restdim; i--){

	vectx[i].rnd(vectx[i]);
	x=vectx[i]; 

	//if (x.sgn() !=0) { 
	if (bounded[i]==0) {
	  
	  lx = x.get_si_exp(expo);

	  nmax=structure[i]+1;
	  
	  // Cf fplll 
	  // Long case 
	  if (expo == 0) {
	    
	    if (lx == 1) {

	      //compteur +=1;
	      
	      somedone = 1;
	      
	      R.subcol(kappa,i,restdim);
	      
	      B.subcol(kappa,i,nmax);
	      
	      if (transf) 
		U.subcol(kappa,i,min(d,nmax));
	      
	    } 
	    else if (lx == -1) {

	      //compteur +=1;
	      
	      somedone = 1;
 
	      R.addcol(kappa,i,restdim);
	      
	      B.addcol(kappa,i,nmax);
	      
	      if (transf) 
		U.addcol(kappa,i,min(d,nmax));

	    } 
	    else { 

	      //compteur +=1;
	      
	      somedone = 1;

	      
	      if (fast_long_flag == 1) {
	      
	       	R.submulcol(kappa,i,x,restdim);
	       	B.addmulcol_si(kappa,i,-lx,nmax);
	       	if (transf)  
	       	  U.addmulcol_si(kappa,i,-lx,min(d,nmax));
		
	      } // end fast_long
	      else {
		
	       	set_f(xz,x);  
	      
	       	R.submulcol(kappa,i,x,restdim);	
	       	B.submulcol(kappa,i,xz,nmax);
	       	if (transf)  
	       	  U.submulcol(kappa,i,xz,min(d,nmax));
	      }	      
	    } 
	    
	  } // end expo == 0 
	  else {  // expo <> 0 

	    //compteur +=1;
	      
	    somedone = 1;
	    
	    // **** À FAIRE Le mettre pour Seysen 
	    if (fast_long_flag == 1) {
	    
	      R.submulcol(kappa,i,x,restdim);
	      B.addmulcol_si_2exp(kappa,i,-lx,expo,nmax);
	      if (transf)  
		U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));
	    
	    } // end fast_long
	    else {
	   
	      set_f(xz,x);  
	  
	      R.submulcol(kappa,i,x,restdim);	
	      B.submulcol(kappa,i,xz,nmax);
	      if (transf)  
		U.submulcol(kappa,i,xz,min(d,nmax));
	  
	    } // end no long


	    
	  } // end expo <> 0 
	} // Non zero combination 

      } // end application de la transformation 

      indexdec+=bdim;     
      ld=ld*2;

    } // End loop on log blocks 

    

    if (somedone) {
       
      compteur+=1;
     
      
      col_kept[kappa]=0;

      t.mul(approx,normB2[kappa]);

      householder_r(kappa);  

      // IICI 
      //nonstop = (normB2[kappa] < t);  // ne baisse quasiment  plus ? 
     
    }
    else {
      
      nonstop=0;
     
     }

  } // end while 

  return somedone;

};




/* -------------------------------------------------------------------------
   Size reduction 

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm

   ------------------------------------------------------------------------- */


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
SLattice<ZT,FT, MatrixZT, MatrixFT>::hsizereduce(int kappa, int fromk) { 

  nmaxkappa=structure[kappa]+1;

  FP_NR<FT> approx;
  
  approx=0.1;


  FP_NR<FT> t,tmpfp;
  Z_NR<ZT>  xz,tmpz;

  long expo,lx;

  int i,w=0;

  bool nonstop=1;
  bool prev_nonstop = 1;
  
  bool somedone=0;

  int nmax; // De la structure triangulaire 
  
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
    
    somedone = 0;

    //compteur +=1;

    // Loop through the column 
    // -----------------------

   
    for (i=startposition; i>-1; i--){  
    
  
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
    
    

    if (somedone) {

      
      compteur+=1;
      
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
	  
      	  FP_NR<FT> theta;
      	  theta = 0.0000001;
      	  theta.mul(theta,R.get(kappa,kappa));
	  
      	  FP_NR<FT> mu,mu_test;

      	  for (i=0; i<kappa; i++) {
	    
      	    mu.div(R.get(i,kappa),R.get(i,i));
      	    mu.abs(mu);

      	    mu_test.div(theta,R.get(i,i));
      	    mu_test.add(mu_test,one);

      	    if (mu.cmp(mu_test) == 1) {
	      
      	      cout << " **** #tests = " << nblov << " **** Anomaly in size reduction, kappa = " << kappa  << endl;
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


/* ------------- */
/* Householder R */
/* ------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
SLattice<ZT,FT, MatrixZT, MatrixFT>::householder_r(int kappa)
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
SLattice<ZT,FT, MatrixZT, MatrixFT>::householder_v(int kappa) 
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
SLattice<ZT,FT, MatrixZT, MatrixFT>::set_nblov_max(unsigned int nb) {

  nblov_max = nb;
  return nblov_max;

}



template<class ZT,class FT,class MatrixZT, class MatrixFT> inline ZZ_mat<ZT>
SLattice<ZT,FT, MatrixZT, MatrixFT>::getbase()
{
   ZZ_mat<ZT> BB(norigin,dorigin);
   for (int i=0; i<norigin; i++) 
     for (int j=0; j<dorigin; j++) BB.Set(i,j,B(i,j)); // reprendre boucle sur les colonnes 

   return BB;
}


  
template<class ZT,class FT, class MatrixZT, class MatrixFT> inline ZZ_mat<ZT>
SLattice<ZT,FT, MatrixZT, MatrixFT>::getU()
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

  
template<class ZT,class FT, class MatrixZT, class MatrixFT> inline  matrix<FP_NR<FT> >
SLattice<ZT,FT, MatrixZT, MatrixFT>::getR()
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
SLattice<ZT,FT, MatrixZT, MatrixFT>::init(int n, int d, bool forU) {

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

  RZ.resize(d,d);
  
  col_kept.resize(d+1); // +1 for the discovery test 
  descendu.resize(d);
  for (i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}

  kappamin.resize(d); // Lowest point down for kappa since last time
  for (j=0; j<d; j++) kappamin[j]=-1;

  VR.resize(d,d);

  U_even.resize(d,d);
  
  U_odd.resize(d,d);

}


template<class ZT,class FT, class MatrixZT, class MatrixFT>
SLattice<ZT,FT, MatrixZT, MatrixFT>::SLattice(ZZ_mat<ZT> A, int S, bool forU, int reduction_method,  int long_flag) {

  // Resizing for dimension divisible by K
  // -------------------------------------
  
  int K=2*S;

  int i,j;

  norigin=A.getRows();
  n=A.getRows();
  dorigin=A.getCols();
  d=A.getCols();

        
  if (d%K !=0) {


    // B.resize(n+K-d%K,d+K-d%K);


    //   for  (i=0; i<K-d%K; i++)
    //     B(i,i)=1;

    //   for  (i=0; i<n; i++)
    //     for (j=0; j<d; j++)
    //       B(i+K-d%K,j+K-d%K)=A(i,j);      

    //   n+=K-d%K;
    //   d+=K-d%K;


      B.resize(n+K-d%K,d+K-d%K);

      Z_NR<mpz_t> tabs,amax;
      amax=0;

      for (i=0; i<n; i++)
      	for (j=0; j<d; j++) {

      	  tabs.abs(A(i,j));
	 
       	  if (tabs.cmp(amax) > 0) amax=tabs;
	 
      	}
      
      for  (i=0; i<n; i++)
        for (j=0; j<d; j++)
          B(i,j)=A(i,j);

      for  (i=0; i<K-d%K; i++)
        B(n+i,d+i)=amax;

      n+=K-d%K;
      d+=K-d%K;

     
    }
  else {
    
    B.resize(n,d);  // Not in init for the mixed matrix case also 

    for (i=0; i<n; i++) 
      for (j=0; j<d; j++) 
	B(i,j)=A.Get(i,j);
    
  }

  // Initializations (after the change of dimension)
  // -----------------------------------------------
  
  init(n,d, forU); 

  if (transf) {    // Not in init for the mixed matrix case also 
    
    U.resize(d,d);
    for (i=0; i<d; i++) U(i,i)=1;     
  }

  nblov_max = 4294967295;
  
  seysen_flag=reduction_method;

  fast_long_flag = long_flag;
  

  matrix_structure(structure, B, n,d);

  for (i=1; i<d; i++)
    structure[i]=max(structure[i-1],structure[i]);
  
  for (i=0; i<d; i++)
    structure[i]=structure[d-1];

  swapstab.resize(K-1);
  for (i=0; i<K-1; i++) 
    swapstab[i]=0;

 }


  
/* --------------------------------------------- */
/* Householder complet */
/* --------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
SLattice<ZT,FT, MatrixZT, MatrixFT>::householder()
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

/* --------------------------------------------- */
/*   Matrix products after block reductions      */
/* --------------------------------------------- */

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
  SLattice<ZT,FT, MatrixZT, MatrixFT>::pmatprod(int S, int dec)
  {
    
    
    int b;
    
    int sdim=d/S; 
    
    // Even case (one could make with only one U)
    // ------------------------------------------
    
    if (dec==0) {

      // Column blocks of B, diag blocks of U_even
      //   simple matrix products 

#ifdef _OPENMP
#pragma omp parallel for 
#endif
      
      for (b=0; b<S; b++) {

	int i,j,k;
	
	int shift = b*sdim;
	  
	MatrixZT tmat;
	tmat.resize(n,sdim);

	for (i=0; i<n; i++) {
	  for (j=shift; j<shift+sdim; j++) {
	    
	    tmat(i,j-shift).mul(B(i,shift),U_even(shift,j));
	    
	    for (k=shift+1; k < shift+sdim; k++) {
	      tmat(i,j-shift).addmul(B(i,k),U_even(k,j));
	    }
	  }
	  for (j=shift; j<shift+sdim; j++)
	    B(i,j)= tmat(i,j-shift);
	} // On rows 
	
	
      } // Block parallel loop 

    }
    // Odd case
    // --------
    else {

      // Column blocks of B, diag blocks of U_even
      //   simple matrix products 

#ifdef _OPENMP
#pragma omp parallel for 
#endif
      
      for (b=0; b<S-1; b++) {
	
	int i,j,k;
	
	int shift = b*sdim+(sdim/2);
	  
	MatrixZT tmat;
	tmat.resize(n,sdim);

	for (i=0; i<n; i++) {
	  for (j=shift; j<shift+sdim; j++) {
	    
	    tmat(i,j-shift).mul(B(i,shift),U_odd(shift,j));
	    
	    for (k=shift+1; k < shift+sdim; k++) {
	      tmat(i,j-shift).addmul(B(i,k),U_odd(k,j));
	    }
	  }
	  for (j=shift; j<shift+sdim; j++)
	    B(i,j)= tmat(i,j-shift);
	} // On rows 
	
	
      } // Block parallel loop 


    } 

      
    
    return 0; 
    
  }

} // end namespace hplll

#endif 

