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
  SLattice<ZT,FT, MatrixZT, MatrixFT>::hlll(double delta, int condbits, int K,  unsigned int lovmax) { 
    
    int bdim;     // Number of blocks and dimension of each block 
                  // Assume that d is a multiple of K >= 4 
                  // K/2 and bdim >= 2 for actual segment) 
    
    bdim = d/K;
    
    int S;   // Number of segments and dimension of each segment  
             // Assume that d is a multiple of K >= 4 
             // K/2 and bdim >= 2 for actual segment) 

    S=K/2;
    
    int Sdim = d/S;

    int i,k;    // block or segment loop 
   
    
        
#ifdef _OPENMP
    OMPTimer time;
    OMPTimer redtime,eventime,oddtime,qrtime,prodtime,esizetime,osizetime,orthotime,totime;
    
    omp_set_num_threads(4);
#else 
    Timer time;
    Timer redtime,eventime,oddtime,qrtime,prodtime,esizetime,osizetime,restsizetime,totime;
#endif 
    
    time.clear();
    redtime.clear();
    eventime.clear();
    oddtime.clear();
    qrtime.clear();
    prodtime.clear();
    esizetime.clear();
    orthotime.clear();
    osizetime.clear();
    totime.clear();
   
    
    int iter;

    bool stop=0;

   
    // *****************************
    // Could be modified dynamically
 
    setprec(condbits);


    // ************************
    // Main loop on block swaps 
    // ************************
   
    totime.start();

    // **** Voir la terminaison avec U Identité
    
    //for (iter=0; iter < 2 ; iter ++){
    for (iter=0; stop==0; iter++) {


      //print2maple(B,n,d);
      
      stop=1;
      
      time.start();

      phouseholder(S);
      
      time.stop();      
      qrtime+=time;
      
   
      // Even block reduction  
      // --------------------

      setId(U);

            
      setId(U_even);
     
      
     
      // The integer block lattice 

      set_f(RZ,R,condbits); 
            
      // DBG
      for (i=0; i<d; i++)
	if (RZ(i,i).sgn() ==0) RZ(i,i)=1;
          
      time.start();

      for (k=0; k<S; k++) { 
  
// #ifdef _OPENMP	
// 	cout << "thread " << omp_get_thread_num() << endl; 
// #endif
	  
	 {
	   cout << "+++++++++++ Even ++++++++++ " << endl;
	   
	  
	   Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  BR(getblock(RZ,k,k,S,0),TRANSFORM,DEF_REDUCTION);
	   // Lattice<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > >  BR(getblock(RZ,k,k,S,0),TRANSFORM,DEF_REDUCTION);

	   BR.set_nblov_max(lovmax);
	   BR.hlll(delta);
	   cout << endl << "even nblov " << BR.nblov << endl; 
	   nblov+=BR.nblov;
	   putblock(RZ,BR.getbase(),k,k,S,0);
	   putblock(U_even,BR.getU(),k,k,S,0);
	}
      }
      
      time.stop();
      redtime+=time; 
      eventime+=time; 

     

// #pragma omp barrier
      
  

	
//       // Size reduction via size reduction of RZ by blocks 
//       // -------------------------------------------------
      
      time.start();

             
       // update blocks above the diagonal 
       even_updateRZ(S);

	
       pmatprod_in(B,U_even,S);

       time.stop();
       prodtime+=time;
       
       stop= (stop &&  isId(U_even));
       //print2maple(getbase(),n,d);
	
       // print2maple(U_even,d,d);
       
       // print2maple(RZ,d,d);


       setId(U_proper);

       // Les blocs diagonaux ne devraient pas avoir à être ré-orthogonalisés

       time.start();
 
       even_hsizereduce(S); // Uproper implicitely updated
       // !!!!! eventuellement que au-dessus de diag pour odd 
       // Both RZ and newRZ are equal 

       time.stop();
       esizetime+=time;
       

       // Si on l'applique à B autant ne pas le faire 3 fois R flottant, RZ et B ?
       
       time.start();

       pmatprod_in(B,U_proper,S);
       

       time.stop();
       prodtime+=time;
             
       stop= (stop &&  isId(U_proper));
       
// //       if (transf) pmatprod_in(Uglob,U,S);
 
       //time.stop();    
       //prodtime+=time; 

       
      // Re-orthogonalization by blocks  !!! should be available from the even lll reductions or the size reduce? 
      // ------------------------------
      // Orthogonalisation partielle diagonale, permet les prochains appels à LLL mais pas un
       //prochain size reduce global 

       time.start();

// #ifdef _OPENMP
// #pragma omp parallel for 
// #endif 
//        for (k=0; k<S-1; k++) {

// 	 int i,j;
	 
// 	 cout << "+++++++++++ Re-ortho ++++++++++ " << endl; 
	 
// 	 {
// 	   ZZ_mat<mpz_t> TR;
// 	   TR.resize(2*Sdim,Sdim+bdim);
// 	   for (i=0; i<2*Sdim; i++) 
// 	     for (j=0; j<Sdim+bdim; j++)
// 	       TR(i,j)=RZ(k*Sdim+i,k*Sdim+j);

	  
// 	   Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(TR);
	   
// 	   for (j=0; j<Sdim+bdim; j++) {
// 	     T.householder_r(j);
// 	     T.householder_v(j);
// 	   } 

// 	   matrix<FP_NR<mpfr_t> > TB;
// 	   TB.resize(Sdim+bdim,Sdim+bdim);

// 	   TB=T.getR();

	 
	  
// 	   matrix<FP_NR<mpfr_t> > TTB;
// 	   TTB.resize(Sdim,Sdim);

// 	   for (i=0; i<Sdim; i++) 
// 	     for (j=0; j<Sdim; j++)
// 	       TTB(i,j)=TB(bdim+i,bdim+j);

// 	   ZZ_mat<mpz_t> TTR;
// 	   TTR.resize(Sdim,Sdim);

	   
// 	   set_f(TTR,TTB,condbits); 
	   
// 	   //DBG
// 	   cout << endl << endl << "**************************************" << d << "  " << Sdim << "   "  << k*Sdim+bdim+Sdim << endl << endl;
	   
// 	   for (i=0; i<Sdim; i++) 
// 	     for (j=0; j<Sdim; j++)
// 	       newRZ(k*Sdim+bdim+i,k*Sdim+bdim+j)=TTR(i,j);

// 	   // DBG
// 	   for (i=0; i<d; i++)
// 	     if (newRZ(i,i).sgn() ==0) newRZ(i,i)=1;
	   
	   
// 	 } // End parallel re-orthogonalization
//        } 


       phouseholder(S);
       set_f(RZ,R,condbits);

       // DBG
       for (i=0; i<d; i++)
	 if (RZ(i,i).sgn() ==0) RZ(i,i)=1;

       time.stop();
       orthotime+=time;
       

	 // Odd block loop 
	 // --------------


       setId(U);
      

       setId(U_odd);

	
       time.start();

    
#ifdef _OPENMP
#pragma omp parallel for 
#endif

      
       
       for (k=0; k<S-1; k++) {

	 cout << "+++++++++++ Odd ++++++++++ " << endl; 
	
	  Lattice<mpz_t, double, matrix<Z_NR<mpz_t> >, matrix<FP_NR<double> > >  BR(getblock(RZ,k,k,S,bdim),TRANSFORM,DEF_REDUCTION);
	  BR.set_nblov_max(lovmax);
	  BR.hlll(delta);
	  cout << endl << "odd nblov " << BR.nblov << endl; 
	  nblov+=BR.nblov;
	  putblock(U_odd,BR.getU(),k,k,S,bdim);	   
	  //putblock(RZ,BR.getbase(),k,k,S,bdim); //Not here: RZ and the orthogonalization were different

	  
 	}
       
       //odd_updateRZ(S); 

       time.stop();
       redtime+=time; 
       oddtime+=time; 

       time.start();
       
       pmatprod_in(B,U_odd,S);
       
       time.stop();
       prodtime+=time;

      stop= (stop &&  isId(U_odd));

      cout << "!!!!!!!!   " << stop << endl;

      
// #pragma omp barrier
      
//       stop=isId(U)*stop;
        
	
//       // Size reduction via size reduction of RZ by blocks 
//       // -------------------------------------------------
      
//       time.start();
      
      pmatprod_in(RZ,U_odd,S);  

//       pmatprod_in(B,U,S);
     

//       if (transf) pmatprod_in(Uglob,U,S);

//       time.stop();    
//       prodtime+=time; 
      
//       // RZ and B same state 
      
//       setId(U);
      
//       time.start();

       setId(U_proper);
       // ICI

       time.start();

       //phouseholder(S);
       

       time.stop();      
       qrtime+=time;
       
       
       //set_f(RZ,R,condbits);
      
       time.start();
      
       odd_hsizereduce(S); // U implicitely updated 

       time.stop();
       osizetime+=time;

       time.start();
       
       pmatprod_in(B,U_proper,S);
       
       time.stop();
       prodtime+=time;
       


//       pmatprod_in(B,U,S);
      
//       if (transf) pmatprod_in(Uglob,U,S);

//       time.stop();    
//       prodtime+=time; 
     
     
      
    } // End main loop: global iterations iter 
    
      totime.stop();

      cout << endl;
      cout << " Householder: " << qrtime << endl;
      cout << " Re-ortho: " << orthotime  << endl;
      cout << " Reductions: " << redtime << endl;
      cout << "   Even reductions: " << eventime << endl;
      cout << "   Odd reductions: " << oddtime << endl;
      cout << " Products: " << prodtime << endl;
      cout << " Even size reds: " << esizetime  << endl;
      cout << " Odd size reds: " << osizetime  << endl;
      cout << " Total time:  " << totime << endl;  
    
    return 0;
    
}


  /* -------------------------------------------------------- */
  /* Update above diagonal after the even phase               */
  /*   since only diagonal blocks have been computed          */
  /* -------------------------------------------------------- */

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void 
  SLattice<ZT,FT, MatrixZT, MatrixFT>::even_updateRZ(int S)  {

    int k;
    
    int Sdim = d/S;

   
      
  // Loop on the block columns 
  // -------------------------
    
  //PPP  
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    // Could use more parallelism 

    for (k=1; k<S ; k++) {

      int i,j,l;
      
      MatrixZT tRZ;
      tRZ.resize(Sdim,Sdim);
      
      MatrixZT  tU;
      tU.resize(Sdim,Sdim);
      
      MatrixZT tres;
      tres.resize(Sdim,Sdim);
      
      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim; j++)
	  tU(i,j)=U_even(k*Sdim+i,k*Sdim+j);
      
      for (l=0; l<k; l++) {
	
	for (i=0; i<Sdim; i++)
	  for (j=0; j<Sdim; j++)
	    tRZ(i,j)=RZ(l*Sdim+i,k*Sdim+j);
	
	matprod(tres,tRZ,tU);
	
	for (i=0; i<Sdim; i++)
	  for (j=0; j<Sdim; j++)
	    RZ(l*Sdim+i,k*Sdim+j)=tres(i,j);
      }
    } 
  }


   /* -------------------------------------------------------- */
  /* Update above diagonal for the odd phase                  */
  /* -------------------------------------------------------- */

  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void 
  SLattice<ZT,FT, MatrixZT, MatrixFT>::odd_updateRZ(int S)  {

    int k;
    
    int Sdim = d/S;

    int bdim = Sdim/2; 
       
  // Loop on the block columns 
  // -------------------------
    
  //PPP  
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    // Could use more parallelism 

    for (k=0; k<S-1 ; k++) {

      int i,j;
      
      MatrixZT  tU;
      tU.resize(Sdim,Sdim);

    
      MatrixZT tRZ;
      tRZ.resize((k+2)*Sdim,Sdim);

      MatrixZT tres;
      tres.resize((k+2)*Sdim,Sdim);
    
      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim; j++)
	  tU(i,j)=U_odd(k*Sdim+bdim+i,k*Sdim+bdim+j);
      
      for (i=0; i<(k+2)*Sdim; i++) 
	for (j=0; j<Sdim; j++)
	  tRZ(i,j)=RZ(i,k*Sdim+bdim+j);
      
	
      matprod(tres,tRZ,tU);

      for (i=0; i<(k+2)*Sdim; i++) 
	for (j=0; j<Sdim; j++)
	  RZ(i,k*Sdim+bdim+j)=tres(i,j);

      
    } 
  }

  
/* -------------------------------------------------------- */
/* Even phase size reduction                                */
/* -------------------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void 
SLattice<ZT,FT, MatrixZT, MatrixFT>::even_hsizereduce(int S)
{


  int k;

  int Sdim = d/S;


  for (int i=0; i<d; i++) 
    for (int j=0; j< d; j++) 
      newRZ(i,j)=RZ(i,j);


  // Loop on the block columns 
  // -------------------------

  //PPP  
#ifdef _OPENMP
#pragma omp parallel for 
#endif 

  for (k=1; k<S ; k++) {

    OMPTimer c;
    c.clear();
    c.start();

    
    int i,j;

   
    // lattice for size reduction 
    ZZ_mat<ZT> tmpM;
    tmpM.resize(Sdim,2*Sdim);
    
    //Lattice<ZT, FT, MatrixZT, MatrixFT> RZloc(tmpM,TRANSFORM,DEF_REDUCTION);
    // Ok since the diagonal blocks are reduced 
    Lattice<ZT, double, MatrixZT, matrix<FP_NR<double> > > RZloc(tmpM,TRANSFORM,DEF_REDUCTION);

    // tmp for U  
    ZZ_mat<ZT> tmpU;
    tmpU.resize(2*Sdim,2*Sdim);

  

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
	  tmpM(i,Sdim+j)=newRZ(l*Sdim+i,k*Sdim+j);


            
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
	  U_proper(l*Sdim+i,k*Sdim+j)=tmpU(i,Sdim+j);

            
      // Update of newRZ for the remaining computations in the block column 
      // Clean matrix product to do  RZ * U in newRZ

            
       for (i=0; i<(l+1)*Sdim; i++)  {  


	int jj,kk;
	
	for (jj=0; jj<Sdim; jj++) {  // i,jj in the result 
	  
	  for (kk=0; kk<Sdim; kk++) {
	    
	    // RZ(i,(l-1)*Sdim+kk)   x  tmpU(kk,Sdim+jj) += newRZ(i,jj) 
	    (newRZ(i,k*Sdim+jj)).addmul( RZ(i,l*Sdim+kk) , tmpU(kk, Sdim+jj) ); 
	  }
	}	
      } // end RZ * U 

             
    } // loop in the block column k 
    
    c.stop();
    
  } // parallel loop on the blocks 


#pragma omp barrier

  // RZ updated at the end since initial value used by threads above
  
  for (int i=0; i<d; i++)
    for (int j=0; j<d; j++) 
      RZ(i,j)=newRZ(i,j);

    


  // Mise a jour de RZ avec newRZ, B et de R à la fin pour le nouveau RZ 

}

/* -------------------------------------------------------- */
/*  Odd phase size reduction                                */
/* -------------------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void 
SLattice<ZT,FT, MatrixZT, MatrixFT>::odd_hsizereduce(int S)
{


  int k;

  int Sdim = d/S;
  int bdim = Sdim/2;

  // Loop on the block columns 
  // -------------------------

  
#ifdef _OPENMP
#pragma omp parallel for 
#endif 


  for (k=1; k<=S ; k++) {

    OMPTimer c;
    c.clear();
    c.start();

    
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

      // ICI
      //cout << "* A ** " << k << "   " << l << endl;
      //print2maple(tmpM,Sdim,2*Sdim-dlast);
		  
      // Local size reduction 

      for (i=0; i<Sdim; i++) {

	RZloc.householder_r(i);
	RZloc.householder_v(i);
      }

      // DBG 
      //print2maple(tmpM,Sdim,2*Sdim-dlast);
 
      for (i=Sdim; i<2*Sdim-dlast; i++) {

	RZloc.hsizereduce(i,Sdim-1);

	
      }

      // Update of U
      tmpU=RZloc.getU();
	
      for (i=0; i<Sdim; i++)
	for (j=0; j<Sdim-dlast; j++) 
	  U_proper(l*Sdim-bdim+i,k*Sdim-bdim+j)=tmpU(i,Sdim+j);


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
    
    tmpM.resize(Sdim,bdim+Sdim-dlast);
    
    Lattice<ZT, FT, MatrixZT, MatrixFT> RZlocf(tmpM,TRANSFORM,DEF_REDUCTION);
    
    // tmp for U  
    tmpU.resize(bdim+Sdim-dlast,bdim+Sdim-dlast);

    // Dummy 
    l=0;

    // Block extraction for size reduction

    for (i=0; i<Sdim; i++)
      for (j=0; j<bdim; j++) 
	tmpM(i,j)=RZ(i,j);

    for (i=0; i<Sdim; i++)
      for (j=0; j<Sdim-dlast; j++) 
	tmpM(i,bdim+j)=newRZ(i,j);
      
    RZlocf.assign(tmpM);

    // ICI
    //cout << "* B ** " << k << "   " << l << endl;
    //print2maple(tmpM,Sdim,bdim+Sdim-dlast);
      
    // Local size reduction 
    
    for (i=0; i<bdim; i++) {
      
      RZlocf.householder_r(i);
      RZlocf.householder_v(i);
    }

    // DBG 
    //print2maple(tmpM,Sdim,bdim+Sdim-dlast);
      
    for (i=bdim; i<bdim+Sdim-dlast; i++) {

      RZlocf.hsizereduce(i,bdim-1);

    }

    // Update of U
    tmpU=RZlocf.getU();
	
    for (i=0; i<bdim; i++)
      for (j=0; j<Sdim-dlast; j++) 
	U_proper(i,k*Sdim-bdim+j)=tmpU(i,bdim+j);

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
SLattice<ZT,FT, MatrixZT, MatrixFT>::phouseholder(int S)
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

      // Zero test 
      if (s.sgn() <=0) 
	V.set(kappa,kappa,1.0); 
      else { 
	V.div(kappa,kappa+1, R.getcol(kappa,kappa+1), s, n-kappa-1);

	nrtmp.div(w,s);
	V.set(kappa,kappa,nrtmp); 
      }

      for(i=kappa+1; i<d; i++)  R.set(i,kappa,0.0); 
       
    }  // end diag computation 

 
    // Parallel application to other blocks 
    // ------------------------------------
    int lb;


#ifdef _OPENMP
#pragma omp parallel for shared (l)
#endif 

    for (lb=l+1; lb<S; lb++) {

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
  ZZ_mat<ZT> BB(norigin,dorigin);
  for (int i=0; i<norigin; i++) 
    for (int j=0; j<dorigin; j++) BB.Set(i,j,B(i,j)); // reprendre boucle sur les colonnes 

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

  newRZ.resize(d,d);
   
   
  U.resize(d,d);
  U_even.resize(d,d);
  U_odd.resize(d,d);
  U_proper.resize(d,d);
 
  if (transf) {
    Uglob.resize(d,d);
    setId(Uglob);
  } 

}


template<class ZT,class FT, class MatrixZT, class MatrixFT>
SLattice<ZT,FT, MatrixZT, MatrixFT>::SLattice(ZZ_mat<ZT> A, int K, bool forU, int reduction_method) {


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
       	  if (tabs.cmp(amax) > 0) amax=A(i,j);
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
  
  init(n,d, forU); 

}




} // end namespace hplll

#endif 

