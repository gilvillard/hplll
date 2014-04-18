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

namespace hplll { 


  // cas rectangles ? 

  template<class ZT,class FT, class MatrixZT, class MatrixFT>  int 
  PLattice<ZT,FT, MatrixZT, MatrixFT>::hlll(double delta) { 

    int K,bdim;   // Number of blocks and dimension of each block 
                  // Assume that d is a multiple of K >= 4 
                  // K/2 and bdim >= 2 for actual segment) 
    K=4;
    bdim = d/K;
    
    int S,sdim;   // Number of segments and dimension of each segment  
                  // Assume that d is a multiple of K >= 4 
                  // K/2 and bdim >= 2 for actual segment) 
    S=K/2;
    sdim = d/S;

    int k;    // block or segment loop 
    int i;

    ZZ_mat<ZT> U;
    U.resize(d,d);

    ZZ_mat<ZT> tmpB;
    tmpB.resize(n,d);
    Lattice<ZT, FT, MatrixZT, MatrixFT> LB(tmpB,TRANSFORM,DEF_REDUCTION);

    int start;
    
    start=utime();
    householder();
    start=utime()-start;
    cout << "   QR " << start/1000 << " ms" << endl;
    
    long condbits;
    

    int iter;

    bool stop=0;

    // Main loop on ...
    // ****************

    for (iter=0; stop==0; iter++) {

      // Even block reduction  
      // --------------------

      setId(U);

      //print2maple(B,n,d);

      condbits=approx_cond();
      cout << endl << "************* Even approx cond " << condbits << endl; 

      set_f(RZ,R,condbits);

      for (k=0; k<S; k++) {   //cout << "+++++++++  Even ++++++++++++ " << endl; 
	 
	//print2maple(getblock(RZ,k,k,S,0),sdim,sdim);

	Lattice<ZT, dpe_t, MatrixZT, MatrixPE<double, dpe_t> > BR(getblock(RZ,k,k,S,0),TRANSFORM,DEF_REDUCTION);
	BR.hlll(delta);
	nblov+=BR.nblov;

	putblock(U,BR.getU(),k,k,S,0);

      }

      stop=isId(U);
      //cout << "Stop: "  <<  stop << endl; 

      matprod_in(B,U);
    
      LB.assign(B);

      for (i=0; i<d; i++) {
	
	LB.hsizereduce(i);
	LB.householder_v(i);
      }

      // Householder have been implicitely computed
      R=LB.getR();
      
      matprod_in(B,LB.getU());


      // Odd block loop 
      // --------------

      setId(U);
      
      condbits=approx_cond();
      cout << endl << "************* Odd approx cond " << condbits << endl;
 
      set_f(RZ,R,condbits);

      for (k=0; k<S-1; k++) {
	//cout << "+++++++++++ Odd ++++++++++ " << endl; 
	
	//print2maple(getblock(RZ,k,k,S,bdim),sdim,sdim);
	
	Lattice<ZT, dpe_t, MatrixZT, MatrixPE<double, dpe_t> > BR(getblock(RZ,k,k,S,bdim),TRANSFORM,DEF_REDUCTION);
	BR.hlll(delta);
	nblov+=BR.nblov;

	putblock(U,BR.getU(),k,k,S,bdim);
      }

      stop=isId(U)*stop;
      //cout << "Stop: "  <<  stop << endl; 
     
      matprod_in(B,U);
    
      LB.assign(B);

      for (i=0; i<d; i++) {
	
	LB.hsizereduce(i);
	LB.householder_v(i);
      }
      
      // Householder have been implicitely computed
      R=LB.getR();
      matprod_in(B,LB.getU());
     
      //print2maple(B,n,d);

    } // End main loop: global iterations iter 




    /*

    // Odd block loop 
    // --------------

    
    setId(U);

    for (k=0; k<S-1; k++) {
      cout << "+++++++++++ Odd ++++++++++ " << endl; 
      
      print2maple(getblock(RZ,k,k,S,bdim),sdim,sdim);
      
      Lattice<ZT, dpe_t, MatrixZT, MatrixPE<double, dpe_t> > BR(getblock(RZ,k,k,S,bdim),TRANSFORM,DEF_REDUCTION);
      BR.hlll(0.99);

      putblock(U,BR.getU(),k,k,S,bdim);
      }

    print2maple(U,d,d);
    
   

    print2maple(LB.getbase(),n,d);*/

  /*start=utime();
  BR.hlll(0.99);
  start=utime()-start;
  cout << "   PLLL " << start/1000 << " ms" << endl;
  
  print2maple(BR.getU(),n/K,n/K);*/

  

  

  return 0;

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
  
  //cout << "**  eta  **  " << eta << endl; 

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

