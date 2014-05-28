/* OpenMP LLL 

Created Mer  9 avr 2014 14:12:40 CEST 
Copyright (C) 2014  Gilles Villard 

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


#ifndef HPLLL_PLLL_CC
#define HPLLL_PLLL_CC

#ifdef _OPENMP
#include <omp.h>
#endif 

namespace hplll { 


  // cas rectangles ? 

  template<class ZT,class FT, class MatrixZT, class MatrixFT>  int 
  PLattice<ZT,FT, MatrixZT, MatrixFT>::hlll(double delta, int K, unsigned int lovmax) { 

    int bdim;     // Number of blocks and dimension of each block 
                  // Assume that d is a multiple of K >= 4 
                  // K/2 and bdim >= 2 for actual segment) 
    
    bdim = d/K;
    
    int S;   // Number of segments and dimension of each segment  
                  // Assume that d is a multiple of K >= 4 
                  // K/2 and bdim >= 2 for actual segment) 
    S=K/2;

    int k;    // block or segment loop 
    int i;


    ZZ_mat<ZT> tmpB;
    tmpB.resize(n,d);
    Lattice<ZT, FT, MatrixZT, MatrixFT> LB(tmpB,NO_TRANSFORM,DEF_REDUCTION);

#ifdef _OPENMP
    OMPTimer time;
    OMPTimer redtime,eventime,oddtime,qrtime,prodtime,esizetime,osizetime,restsizetime;

    omp_set_num_threads(4);
#else 
    Timer time;
    Timer redtime,eventime,oddtime,qrtime,prodtime,esizetime,osizetime,restsizetime;
#endif 
    
    time.clear();
    redtime.clear();
    eventime.clear();
    oddtime.clear();
    qrtime.clear();
    prodtime.clear();
    esizetime.clear();
    restsizetime.clear();
    osizetime.clear();

    time.start();

    householder();

    time.stop();
    qrtime+=time; 



    long condbits;
    
    int iter;

    bool stop=0;



    // ************************
    // Main loop on block swaps 
    // ************************
    // Size reduced in input 
    

    for (iter=0; stop==0; iter++) {
      //for (iter=0; iter < 1 ; iter ++){

      // Even block reduction  
      // --------------------

      setId(U);

      condbits=approx_cond();
      cout << endl << "************* Even approx cond " << condbits << "    " << "S = " << S << endl; 
      cout << " Reductions: " << redtime << endl;
      cout << " Even reductions:   " << eventime << endl;
      cout << " Odd reductions:   " << oddtime << endl;
      cout << " Products:   " << prodtime << endl;
      cout << " Calls even size reds:  " << esizetime << endl;
      cout << " Rest even size reds:  " << restsizetime << endl;

      LB.setprec(condbits);

      set_f(RZ,R,condbits);  // Le limiter aux blocs 
      
      time.start();

#ifdef _OPENMP
#pragma omp parallel for 
#endif 
      
      for (k=0; k<S; k++) {   
#ifdef _OPENMP	
	cout << "thread " << omp_get_thread_num() << endl; 
#endif
	// Double hence no global prec problem 
	Lattice<ZT, dpe_t, MatrixZT, MatrixPE<double, dpe_t> > BR(getblock(RZ,k,k,S,0),TRANSFORM,DEF_REDUCTION);
	BR.set_nblov_max(lovmax);
	BR.hlll(delta);
	cout << endl << "even nblov " << BR.nblov << endl; 
	nblov+=BR.nblov;
	putblock(U,BR.getU(),k,k,S,0);

      }
      
#ifdef _OPENMP
#pragma omp barrier
#endif 

      
      time.stop();
      redtime+=time; 
      eventime+=time; 

      stop=isId(U);

      int method=0;

      if (method == 1) {

	// Initial even size reduction method via a lattice with via size reduction of RZ by blocks 
	// -------------------------------------------------
	
	time.start();

	matprod_in(B,U);

	time.stop();    
	prodtime+=time; 

	LB.assign(B);

	time.start();

	for (i=0; i<d; i++) {
	
	  LB.hsizereduce(i);
	  LB.householder_v(i);
	}

	time.stop();
	esizetime+=time;


	R=LB.getR();
      
	set(B,LB.getbase());

      } // end method 0 for size reduction 
      else { // new method for even size reduction 
	
	// Size reduction via size reduction of RZ by blocks 
	// -------------------------------------------------

	matprod_in(RZ,U);  

	time.start();

	matprod_in(B,U);

	time.stop();    
	prodtime+=time; 

	// RZ and B same state 

	setId(U);

	time.start();

	even_hsizereduce(S,condbits); // U implicitely updated 

	time.stop();
	esizetime+=time;

	



	matprod_in(B,U);

	time.start();

	// CHANGER LE PREC MPFR EN FONCTION DE RZ CAR A BAISSÉE
	LB.assign(B);
	for (i=0; i<d; i++) {
	  
	  LB.householder_r(i);
	  LB.householder_v(i);
	}

	

	time.stop();
	restsizetime+=time;

	R=LB.getR();

	// INUTILE SANS DOUTE 
	set(B,LB.getbase());

      } // end new method for even size reduction 

      // Odd block loop 
      // --------------

      setId(U);
      
      condbits=approx_cond();
    
      set_f(RZ,R,condbits);

      time.start();

#ifdef _OPENMP
#pragma omp parallel for 
#endif 
      for (k=0; k<S-1; k++) {
	//cout << "+++++++++++ Odd ++++++++++ " << endl; 
	
	Lattice<ZT, dpe_t, MatrixZT, MatrixPE<double, dpe_t> > BR(getblock(RZ,k,k,S,bdim),TRANSFORM,DEF_REDUCTION);
	BR.set_nblov_max(lovmax);
	BR.hlll(delta);
	cout << endl << "odd nblov " << BR.nblov << endl;
	nblov+=BR.nblov;

	putblock(U,BR.getU(),k,k,S,bdim);
      }

      
      time.stop();
      redtime += time;
      oddtime += time;

      stop=isId(U)*stop;
           
      time.start();

      matprod_in(B,U);

      time.stop();
      prodtime+=time;


      LB.assign(B);


      time.start();

      for (i=0; i<d; i++) {
	
	LB.hsizereduce(i);
	LB.householder_v(i);
      }
      
      time.stop();
      osizetime+=time;


      // Householder have been implicitely computed

      R=LB.getR();  
      set(B,LB.getbase());


    } // End main loop: global iterations iter 

    cout << endl;
    cout << " Initial QR  " << qrtime << endl;
    cout << " Reductions: " << redtime << endl;
    cout << " Even reductions:   " << eventime << endl;
    cout << " Odd reductions:   " << oddtime << endl;
    cout << " Products:   " << prodtime << endl;
    cout << " Calls even size reds:  " << esizetime  << endl;
    cout << " Rest even size reds:  " << restsizetime  << endl;
    cout << " Odd size reds:  " << osizetime  << endl;
  return 0;

}

/* -------------------------------------------------------- */
/* Even phase size reduction                                */
/* -------------------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void 
PLattice<ZT,FT, MatrixZT, MatrixFT>::even_hsizereduce(int S, int prec)
{


  int k;

  int Sdim = d/S;

  // Loop on the block columns 
  // -------------------------

  
#ifdef _OPENMP
#pragma omp parallel for shared (prec)
#endif 

  for (k=1; k<S ; k++) {

    OMPTimer c;
    c.clear();
    c.start();

    mpfr_set_default_prec(prec);

    int i,j;

   
    // lattice for size reduction 
    ZZ_mat<ZT> tmpM;
    tmpM.resize(Sdim,2*Sdim);
    
    Lattice<ZT, FT, MatrixZT, MatrixFT> RZloc(tmpM,TRANSFORM,DEF_REDUCTION);
    // Ok since the diagonal blocks are reduced 
    //Lattice<ZT, dpe_t, MatrixZT, MatrixPE<double, dpe_t> > RZloc(tmpM,TRANSFORM,DEF_REDUCTION);


    // tmp for U  
    ZZ_mat<ZT> tmpU;
    tmpU.resize(2*Sdim,2*Sdim);

    // The local slice of RZ and update 
    MatrixZT newRZ;
    newRZ.resize(d,Sdim);

    for (i=0; i<d; i++) 
      for (j=0; j< Sdim; j++) 
	newRZ(i,j)=RZ(i,k*Sdim+j);

    
    // Loop in the block column 
    // ------------------------

    int l;

    for (l=k-1; l>-1; l--) {

      
      // Block extraction for size reduction
      // from RZ and the update in newR 

      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim; j++) 
	  tmpM(i,j)=RZ(l*Sdim+i,l*Sdim+j);

      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim; j++) 
	  tmpM(i,Sdim+j)=newRZ(l*Sdim+i,j);

      
      RZloc.assign(tmpM);
 
     
     
      // Local size reduction 

      for (i=0; i<Sdim; i++) {

	RZloc.householder_r(i);
	RZloc.householder_v(i);
      }


      for (i=Sdim; i<2*Sdim; i++) {

	RZloc.hsizereduce(i,Sdim-1);

	
      }

      // Update of U
      tmpU=RZloc.getU();

      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim; j++) 
	  U(l*Sdim+i,k*Sdim+j)=tmpU(i,Sdim+j);

      // Update of newRZ for the remaining computations in the block column 
      // Clean matrix product to do  RZ * U in newRZ

      for (i=0; i<l*Sdim; i++)  {  
 
	int jj,kk;
	Z_NR<ZT> p;
	
	for (jj=0; jj<Sdim; jj++) {  // i,jj in the result 
	  p =0; 
	  for (kk=0; kk<Sdim; kk++) {
	    
	    // RZ(i,(l-1)*Sdim+kk)   x  tmpU(kk,Sdim+jj) += newRZ(i,jj) 
	    (newRZ(i,jj)).addmul( RZ(i,l*Sdim+kk) , tmpU(kk, Sdim+jj) ); 
	  }
	}	
      } // end RZ * U 

      
    } // loop in the block column 

    c.stop();
    cout << "On block : " << k << "  " << c << endl; 
  } // parallel loop on the blocks 



  // Mise a jour de RZ avec newRZ, B et de R à la fin pour le nouveau RZ 

}

/* -------------------------------------------------------- */
/* Approximate log_2 of condition number of R, size-reduced */
/* -------------------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline long 
PLattice<ZT,FT, MatrixZT, MatrixFT>::approx_cond()
{

  int i,j;

  FP_NR<FT>  tmp,ttmp;

  // Size-reducedness 
  // ----------------
  
  FP_NR<FT>  eta;
  eta=0.0;

  for  (j=1; j<d; j++) 

    for (i=0; i<j; i++) {
      tmp.div(R(i,j),R(i,i));
      tmp.abs(tmp);
      if (tmp.cmp(eta) > 0) eta=tmp;
    }
  
  // cout << "**  eta  **  " << eta << endl; 

  // diag quo
  // --------
  
  FP_NR<FT>  maxquo;
  maxquo=0.0;

  for  (j=0; j<d; j++) 
    for (i=0; i<d; i++) {
      tmp.div(R(j,j),R(i,i));
      tmp.abs(tmp);
      if (tmp.cmp(maxquo) > 0) maxquo=tmp;
    }
  
  //cout << "**  maxquo  **  " << maxquo << endl;

  // Approx cond all together 
  // ------------------------

  long condb; 
 
  condb=maxquo.exponent();
  tmp=1.0;
  eta.add(tmp,eta);
  condb += (d-1) + eta.exponent();

  Z_NR<long> dd;
  dd=d;
  condb += 1 + dd.exponent();

  return condb;
}



/* --------------------------------------------- */
/* Complete Householder */
/* --------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
PLattice<ZT,FT, MatrixZT, MatrixFT>::householder()
{

  int i,k,kappa;
  FP_NR<FT> nrtmp,s,w; 
  

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





template<class ZT,class FT,class MatrixZT, class MatrixFT> inline ZZ_mat<ZT> PLattice<ZT,FT, MatrixZT, MatrixFT>::getbase()
{
  ZZ_mat<ZT> BB(n,d);
  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) BB.Set(i,j,B(i,j)); // reprendre boucle sur les colonnes 

  return BB;
}



template<class ZT,class FT, class MatrixZT, class MatrixFT> inline  matrix<FP_NR<FT> > PLattice<ZT,FT, MatrixZT, MatrixFT>::getR()
{
  matrix<FP_NR<FT> >  RR(d,d);
  FP_NR<FT> tmp;

  for (int i=0; i<d; i++) 
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
PLattice<ZT,FT, MatrixZT, MatrixFT>::init(int n, int d) {

  nblov=0;

  R.resize(n,d);

  V.resize(n,d);

  RZ.resize(d,d);

  U.resize(d,d);

}


template<class ZT,class FT, class MatrixZT, class MatrixFT>
PLattice<ZT,FT, MatrixZT, MatrixFT>::PLattice(ZZ_mat<ZT> A) {

  
  n=A.getRows();
  d=A.getCols();

  init(n,d); 

  int i,j;

  B.resize(n,d);  // Not in init for the mixed matrix case also 

  for (i=0; i<n; i++) 
    for (j=0; j<d; j++) 
      B(i,j)=A.Get(i,j);

 }




} // end namespace hplll

#endif 

