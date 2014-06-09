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


#ifndef HPLLL_SLLL_CC
#define HPLLL_SLLL_CC

#ifdef _OPENMP
#include <omp.h>
#endif 

namespace hplll { 


  // cas rectangles ? 

  template<class ZT,class FT, class MatrixZT, class MatrixFT>  int 
  SLattice<ZT,FT, MatrixZT, MatrixFT>::hlll(double delta, int K, int level, unsigned int lovmax) { 
    
    int bdim;     // Number of blocks and dimension of each block 
                  // Assume that d is a multiple of K >= 4 
                  // K/2 and bdim >= 2 for actual segment) 
    
    bdim = d/K;
    
    int S;   // Number of segments and dimension of each segment  
             // Assume that d is a multiple of K >= 4 
             // K/2 and bdim >= 2 for actual segment) 

    int SP=K/2;

    S=K/2;
    
    int k;    // block or segment loop 
   
    
        
#ifdef _OPENMP
    OMPTimer time;
    OMPTimer redtime,eventime,oddtime,qrtime,prodtime,esizetime,osizetime,restsizetime,totime,ttime;
    
    omp_set_num_threads(4);
#else 
    Timer time;
    Timer redtime,eventime,oddtime,qrtime,prodtime,esizetime,osizetime,restsizetime,totime,ttime;
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
    totime.clear();
    ttime.clear();
   

    long condbits;
    
    int iter;

    bool stop=0;


    // The input matrix is assumed to be triangular and size reduced 
    // -------------------------------------------------------------
    
    Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > >  LI(getbase());
    condbits = LI.lcond(TRIANGULAR_PROPER);
    
    setprec(condbits);

    time.start();

    householder();

    time.stop();
    cout << "++++++++++++ Level: " << level << "  +++++++++++++++++++ Prec: " << mpfr_get_default_prec() << "   " << time << endl;  
    qrtime+=time; 



    
    // ************************
    // Main loop on block swaps 
    // ************************
    // Size reduced in input 
    

    totime.start();

    for (iter=0; stop==0; iter++) {
      //          for (iter=0; iter < 1 ; iter ++){

      
      time.start();

      phouseholder(SP,condbits);
            
      time.stop();      
      restsizetime+=time;

      

      // Even block reduction  
      // --------------------

      setId(U);

      condbits=approx_cond();

      cout << endl << "******* Level: " << level << "  ****** Even approx cond " << condbits << "    " << "S = " << S << endl; 
          
      cout << " Reductions: " << redtime << endl;
      //cout << " Even reductions:   " << eventime << endl;
      //cout << " Odd reductions:   " << oddtime << endl;
      cout << " Products:   " << prodtime << endl;
      cout << " Calls even size reds:  " << esizetime << endl;
      cout << " Rest size reds:  " << restsizetime << endl;
      totime.stop();
      ttime+=totime;
      cout << " Total time:  " << ttime << endl;  
      cout << " Nblov : " << nblov << endl; 
      totime.start();
       
     

      set_f(RZ,R,condbits);  // Le limiter aux blocs 
      
      setprec(condbits);
   
      time.start();

#ifdef _OPENMP
#pragma omp parallel for shared(condbits)
#endif 
      
      for (k=0; k<S; k++) {   
#ifdef _OPENMP	
	//cout << "thread " << omp_get_thread_num() << endl; 
#endif
	// Double hence no global prec problem 
 
	if (level ==0) {
	  Lattice<ZT, dpe_t, MatrixZT, MatrixPE<double, dpe_t> > BR(getblock(RZ,k,k,S,0),TRANSFORM,DEF_REDUCTION);
	  BR.set_nblov_max(lovmax);
	  BR.hlll(delta);
	  cout << endl << "even nblov " << BR.nblov << endl; 
	  nblov+=BR.nblov;
	  putblock(U,BR.getU(),k,k,S,0);
	} 
	else {
	  mpfr_set_default_prec(condbits);
	  SLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > BP(getblock(RZ,k,k,S,0),TRANSFORM,DEF_REDUCTION);
	  BP.hlll(delta,4,level-1);
	  cout << " ========= ICI " << endl; 
	  //cout << endl << "even nblov " << BP.nblov << endl; 
	  nblov+=BP.nblov;
	  putblock(U,BP.getU(),k,k,S,0);
	  mpfr_set_default_prec(condbits); 
	}
      }
      
      mpfr_set_default_prec(condbits); 

      time.stop();
      redtime+=time; 
      eventime+=time; 
      
      stop=isId(U);

	
      // Size reduction via size reduction of RZ by blocks 
      // -------------------------------------------------
      
      time.start();
      
      pmatprod_in(RZ,U,SP);  
      
      pmatprod_in(B,U,SP);
      

      if (transf) pmatprod_in(Uglob,U,SP);
 
      time.stop();    
      prodtime+=time; 
      
      // RZ and B same state 
      
      setId(U);
      
      time.start();
      
      bool refresh = true; // False: computes Rt 
 

      even_hsizereduce(S,condbits,refresh); // U implicitely updated 
      
      time.stop();
      esizetime+=time;
	
      time.start();
      
      pmatprod_in(B,U,SP);


      if (transf) pmatprod_in(Uglob,U,SP);
 
      time.stop();    
      prodtime+=time; 
      
      time.start();
      
      // CHANGER LE PREC MPFR EN FONCTION DE RZ CAR A BAISSÉE
      if (refresh) 
	phouseholder(SP,condbits);
      //else  
      //set(R,Rt);
      

      time.stop();
      //cout << "+++++++++++++++++++++++++++++++ Prec: " << mpfr_get_default_prec() << "   " << time << endl;
      restsizetime+=time;
      
    
      // Odd block loop 
      // --------------
      
      setId(U);
      
      condbits=approx_cond();

      
      set_f(RZ,R,condbits);
     

      setprec(condbits);
  

      time.start();

#ifdef _OPENMP
#pragma omp parallel for 
#endif 
      for (k=0; k<S-1; k++) {
	//cout << "+++++++++++ Odd ++++++++++ " << endl; 
	
	if (level ==0) {
	  Lattice<ZT, dpe_t, MatrixZT, MatrixPE<double, dpe_t> > BR(getblock(RZ,k,k,S,bdim),TRANSFORM,DEF_REDUCTION);
	  //BR.set_nblov_max(lovmax);
	  BR.hlll(delta);
	  //cout << endl << "odd nblov " << BR.nblov << endl;
	  nblov+=BR.nblov;

	  putblock(U,BR.getU(),k,k,S,bdim);
	} 
	else {
	  mpfr_set_default_prec(condbits);
	  SLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > BP(getblock(RZ,k,k,S,bdim),TRANSFORM,DEF_REDUCTION);
	  BP.hlll(delta,4,level-1);
	  cout << " ========= ICI " << endl; 
	  //cout << endl << "even nblov " << BP.nblov << endl; 
	  nblov+=BP.nblov;
	  putblock(U,BP.getU(),k,k,S,bdim);
	  mpfr_set_default_prec(condbits); 
	}
      }
      mpfr_set_default_prec(condbits); 
      
      time.stop();
      redtime += time;
      oddtime += time;

      stop=isId(U)*stop;
        
	
      // Size reduction via size reduction of RZ by blocks 
      // -------------------------------------------------
      
      time.start();
      
      pmatprod_in(RZ,U,SP);  

      pmatprod_in(B,U,SP);
     

      if (transf) pmatprod_in(Uglob,U,SP);

      time.stop();    
      prodtime+=time; 
      
      // RZ and B same state 
      
      setId(U);
      
      time.start();
      
      odd_hsizereduce(S,condbits); // U implicitely updated 
      
      time.stop();
      osizetime+=time;
      
      time.start();


      pmatprod_in(B,U,SP);
      
      if (transf) pmatprod_in(Uglob,U,SP);

      time.stop();    
      prodtime+=time; 
     
     
      
    } // End main loop: global iterations iter 
    
    cout << endl;
    cout << " Initial QR  " << qrtime << endl;
    cout << " Reductions: " << redtime << endl;
    cout << " Even reductions:   " << eventime << endl;
    cout << " Odd reductions:   " << oddtime << endl;
    cout << " Products:   " << prodtime << endl;
    cout << " Calls even size reds:  " << esizetime  << endl;
    cout << " Calls odd size reds:  " << osizetime  << endl;
    cout << " Rest size reds:  " << restsizetime  << endl;
    totime.stop();
    ttime+=totime;
    cout << " Total time:  " << ttime << endl;  

    totime.start();
    
    return 0;
    
}
  
/* -------------------------------------------------------- */
/* Even phase size reduction                                */
/* -------------------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void 
SLattice<ZT,FT, MatrixZT, MatrixFT>::even_hsizereduce(int S, int prec, bool refresh)
{


  int k;

  int Sdim = d/S;

  // Loop on the block columns 
  // -------------------------

  //PPP  
#ifdef _OPENMP
#pragma omp parallel for shared (prec,refresh)
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

    matrix<FP_NR<FT> > tmpR;
    tmpR.resize(Sdim,2*Sdim);
    
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

      if (refresh == false) { // Householder need to be available via RZ 

	tmpR=RZloc.getR();
	for (i=0; i<Sdim; i++)
	  for (j=0; j<Sdim; j++) {
	  
	    Rt.set(l*Sdim+i,l*Sdim+j,tmpR(i,j)); 
	  Rt.set(l*Sdim+i,k*Sdim+j,tmpR(i,Sdim+j));
	} 
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
	
	for (jj=0; jj<Sdim; jj++) {  // i,jj in the result 
	  
	  for (kk=0; kk<Sdim; kk++) {
	    
	    // RZ(i,(l-1)*Sdim+kk)   x  tmpU(kk,Sdim+jj) += newRZ(i,jj) 
	    (newRZ(i,jj)).addmul( RZ(i,l*Sdim+kk) , tmpU(kk, Sdim+jj) ); 
	  }
	}	
      } // end RZ * U 

      
    } // loop in the block column 


    // Last diagonal block of Householder may be required 
    if ((refresh == false) && (k==(S-1))) { 

      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim; j++) 
	  tmpM(i,j)=RZ(k*Sdim+i,k*Sdim+j);
      
      
      RZloc.assign(tmpM);
 
      // Local size reduction 

      for (i=0; i<Sdim; i++) {

	RZloc.householder_r(i);
	RZloc.householder_v(i);
      }

      // 
      tmpR=RZloc.getR();
      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim; j++) {
	  
	  Rt.set(k*Sdim+i,k*Sdim+j,tmpR(i,j)); 
	 
	} 
    } // Last block for householder if not refreshed 
    
    
    c.stop();
    //cout << "On block : " << k << "  " << c << endl; 
  } // parallel loop on the blocks 



  // Mise a jour de RZ avec newRZ, B et de R à la fin pour le nouveau RZ 

}

/* -------------------------------------------------------- */
/*  Odd phase size reduction                                */
/* -------------------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void 
SLattice<ZT,FT, MatrixZT, MatrixFT>::odd_hsizereduce(int S, int prec)
{


  int k;

  int Sdim = d/S;
  int bdim = Sdim/2;

  // Loop on the block columns 
  // -------------------------

  //PPP  
#ifdef _OPENMP
  #pragma omp parallel for shared (prec)
#endif 


  for (k=1; k<=S ; k++) {

    OMPTimer c;
    c.clear();
    c.start();

    mpfr_set_default_prec(prec);

    int i,j;

    int dlast=0;

    // Sdim blocks and last smaller block 
    // ----------------------------------

    if (k==S) 
      dlast = bdim;
    else 
      dlast = 0;
 
    // lattice for size reduction 
    ZZ_mat<ZT> tmpM;
    tmpM.resize(Sdim,2*Sdim-dlast);
    
    Lattice<ZT, FT, MatrixZT, MatrixFT> RZloc(tmpM,TRANSFORM,DEF_REDUCTION);
      
    // tmp for U  
    ZZ_mat<ZT> tmpU;
    tmpU.resize(2*Sdim-dlast,2*Sdim-dlast);

    // The local slice of RZ and update 
    MatrixZT newRZ;
    newRZ.resize(d,Sdim-dlast);

    for (i=0; i<d; i++) 
      for (j=0; j< Sdim-dlast; j++) 
	newRZ(i,j)=RZ(i,k*Sdim-bdim+j);

    
    // Loop in the block column 
    // ------------------------

    int l;

    for (l=k-1; l>0; l--) { 
      
      // Block extraction for size reduction
      // from RZ and the update in newR 
      
      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim; j++) 
	  tmpM(i,j)=RZ(l*Sdim-bdim+i,l*Sdim-bdim+j);
      
      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim-dlast; j++) 
	  tmpM(i,Sdim+j)=newRZ(l*Sdim-bdim+i,j);
      
      
      RZloc.assign(tmpM);
 
      // Local size reduction 

      for (i=0; i<Sdim; i++) {

	RZloc.householder_r(i);
	RZloc.householder_v(i);
      }


      for (i=Sdim; i<2*Sdim-dlast; i++) {

	RZloc.hsizereduce(i,Sdim-1);

	
      }

      // Update of U
      tmpU=RZloc.getU();
	
      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim-dlast; j++) 
	  U(l*Sdim-bdim+i,k*Sdim-bdim+j)=tmpU(i,Sdim+j);


      // Update of newRZ for the remaining computations in the block column 
      // Clean matrix product to do  RZ * U in newRZ

      for (i=0; i<l*Sdim-bdim; i++)  {  
 
	int jj,kk;
  	
	for (jj=0; jj<Sdim-dlast; jj++) {  // i,jj in the result 
	   
	  for (kk=0; kk<Sdim; kk++) {
	    
	    // RZ(i,(l-1)*Sdim+kk)   x  tmpU(kk,Sdim+jj) += newRZ(i,jj) 
	    (newRZ(i,jj)).addmul( RZ(i,l*Sdim-bdim+kk) , tmpU(kk, Sdim+jj) ); 
	  }
	}	
      } // end RZ * U 
     
    } // end loop in the block column 
    
    // With respect to the first smaller block l=0
    // -------------------------------------------
    
    // lattice for size reduction 
    
    tmpM.resize(bdim,bdim+Sdim-dlast);
    
    Lattice<ZT, FT, MatrixZT, MatrixFT> RZlocf(tmpM,TRANSFORM,DEF_REDUCTION);
    
    // tmp for U  
    tmpU.resize(bdim+Sdim-dlast,bdim+Sdim-dlast);

    // Dummy 
    l=0;

    // Block extraction for size reduction

    for (i=0; i<bdim; i++)
      for (j=0; j<bdim; j++) 
	tmpM(i,j)=RZ(i,j);

    for (i=0; i<bdim; i++)
      for (j=0; j<Sdim-dlast; j++) 
	tmpM(i,bdim+j)=newRZ(i,j);
      
    RZlocf.assign(tmpM);
 
    // Local size reduction 
    
    for (i=0; i<bdim; i++) {
      
      RZlocf.householder_r(i);
      RZlocf.householder_v(i);
    }
    
    for (i=bdim; i<bdim+Sdim-dlast; i++) {

      RZlocf.hsizereduce(i,bdim-1);

    }

    // Update of U
    tmpU=RZlocf.getU();
	
    for (i=0; i<bdim; i++)
      for (j=0; j<Sdim-dlast; j++) 
	U(i,k*Sdim-bdim+j)=tmpU(i,bdim+j);

    // End w.r.t the first block 
    
    
    c.stop();
    //cout << "On block : " << k << "  " << c << endl; 
    
  } // end parallel loop on the blocks 

}

/* -------------------------------------------------------- */
/* Approximate log_2 of condition number of R, size-reduced */
/* -------------------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline long 
SLattice<ZT,FT, MatrixZT, MatrixFT>::approx_cond()
{

  int i,j;

  FP_NR<FT>  tmp,ttmp;

  // Size-reducedness 
  // ----------------
  
  FP_NR<FT>  eta;
  eta=0.0;

  for  (j=1; j<d; j++) 

    for (i=0; i<j; i++) {
      tmp.div(R.get(i,j),R.get(i,i));
      tmp.abs(tmp);
      if (tmp.cmp(eta) > 0) eta=tmp;
    }
  
  //cout << "************************************************************  eta  **  " << eta << endl; 

  // diag quo
  // --------
  
  FP_NR<FT>  maxquo;
  maxquo=0.0;

  for  (j=0; j<d; j++) 
    for (i=0; i<d; i++) {
      tmp.div(R.get(j,j),R.get(i,i));
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
SLattice<ZT,FT, MatrixZT, MatrixFT>::householder()
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

/* --------------------------------------------- */
/* Complete Householder */
/* --------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
SLattice<ZT,FT, MatrixZT, MatrixFT>::phouseholder(int S, int prec)
{

  int i,k,kappa;
  FP_NR<FT> nrtmp,s,w; 

 
  int Sdim=d/S;

  int l;

  for (l=0; l<S; l++)   {

    // Diagonal block
    // --------------

    for (kappa=l*Sdim; kappa<(l+1)*Sdim; kappa++) {

      if (l==0) 
	R.setcol(kappa,B.getcol(kappa),0,n);

      for (k=l*Sdim; k<kappa; k++) {
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
       
    }  // end diag computation 

 
    // Parallel application to other blocks 
    // ------------------------------------
    int lb;


#ifdef _OPENMP
#pragma omp parallel for shared (prec,l)
#endif 

    for (lb=l+1; lb<S; lb++) {

      mpfr_set_default_prec(prec);

      int kk; 
      int pkappa; 
      FP_NR<FT> pnrtmp;
      
      for (pkappa=lb*Sdim; pkappa<(lb+1)*Sdim; pkappa++) {

	if (l==0) 
	  R.setcol(pkappa,B.getcol(pkappa),0,n);
      		
	for (kk=l*Sdim; kk<(l+1)*Sdim; kk++) {
	  scalarprod(pnrtmp, V.getcol(kk,kk), R.getcol(pkappa,kk), n-kk);
	  R.fmasub(pkappa,kk,R.getcol(pkappa,kk), V.getcol(kk,kk), pnrtmp, n-kk); 
	}
      }
    } 


  } // end block loop 
   
    return 0; 
}






template<class ZT,class FT,class MatrixZT, class MatrixFT> inline ZZ_mat<ZT> SLattice<ZT,FT, MatrixZT, MatrixFT>::getbase()
{
  ZZ_mat<ZT> BB(n,d);
  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) BB.Set(i,j,B(i,j)); // reprendre boucle sur les colonnes 

  return BB;
}



template<class ZT,class FT, class MatrixZT, class MatrixFT> inline  matrix<FP_NR<FT> > SLattice<ZT,FT, MatrixZT, MatrixFT>::getR()
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


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline ZZ_mat<ZT> SLattice<ZT,FT, MatrixZT, MatrixFT>::getU()
{
 
  if (transf) { 
    ZZ_mat<ZT> UU(d,d); 
    for (int i=0; i<d; i++) 
      for (int j=0; j<d; j++) UU.Set(i,j,Uglob.get(i,j)); // reprendre boucle sur les colonnes 
    return UU;
  }
  else {
    cout << "*** Error, PLLL, the transformation matrix has not been computed" << endl;
    ZZ_mat<ZT> UU(0,0);
    return UU;
  }
}

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline unsigned int 
SLattice<ZT,FT, MatrixZT, MatrixFT>::getprec() {

  return (R.get(0,0)).getprec(); 
  
}


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline unsigned int 
SLattice<ZT,FT, MatrixZT, MatrixFT>::setprec(unsigned int prec) {

  // Re-initializations
  // Should be for each variable, non global? 
  // push and pop ? 
 
  unsigned oldprec;
  oldprec=getprec();

  mpfr_set_default_prec(prec);

  R.clear();
  R.resize(n,d);

  Rt.clear();
  Rt.resize(n,d);

  unsigned newprec;
  newprec=(R.get(0,0)).getprec(); 
  if (newprec == oldprec) cout << "Warning: in function setprec plll, the change of precision has no effect" << endl; 

  V.clear();
  V.resize(n,d);

  // The old precision 
  return oldprec;
}


// Constructeur 
// ------------


template<class ZT,class FT, class MatrixZT, class MatrixFT> void 
SLattice<ZT,FT, MatrixZT, MatrixFT>::init(int n, int d, bool forU) {

  transf = forU;

  nblov=0;

  R.resize(n,d);

  Rt.resize(n,d);

  V.resize(n,d);

  RZ.resize(d,d);

  U.resize(d,d);

  if (transf) {
    Uglob.resize(d,d);
    setId(Uglob);
  } 

}


template<class ZT,class FT, class MatrixZT, class MatrixFT>
SLattice<ZT,FT, MatrixZT, MatrixFT>::SLattice(ZZ_mat<ZT> A, bool forU, int reduction_method) {

  

  n=A.getRows();
  d=A.getCols();

  init(n,d, forU); 

  int i,j;

  B.resize(n,d);  // Not in init for the mixed matrix case also 

  for (i=0; i<n; i++) 
    for (j=0; j<d; j++) 
      B(i,j)=A.Get(i,j);

 }




} // end namespace hplll

#endif 

