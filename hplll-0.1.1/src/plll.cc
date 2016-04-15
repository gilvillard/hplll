/* Householder LLL 

Created Lun  9 jui 2014 15:04:23 CEST  
Copyright (C) 2014      Gilles Villard 

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


// MPFR setprec : isreduced and cond are re-initializing the floating point matrices 

namespace hplll { 



  template<class ZT,class FT, class MatrixZT, class MatrixFT>  int 
  PLattice<ZT,FT, MatrixZT, MatrixFT>::hlll(double delta, bool verbose) { 

  
    int kappa=0,i,j,K;
  
    vector<FP_NR<FT> >  prevR(d);


    FP_NR<FT> newt; //testaccu;

    FP_NR<FT> deltab,lovtest;
    deltab=delta;   // TO SEE 

    FP_NR<FT> tmpswap;

    FP_NR<FT> s,sn; // newt test 

    for (i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}


    int odd=0;

    int swapin=1;

    Timer time, tsize,tup;
    time.clear();
    tsize.clear();
    tup.clear();
  
    ichrono.clear();
    itime.clear();
  
  
    // Main iterative loop 
    // -------------------

    int endposition=d; 
      
    for (K=0; (swapin > 0) || (odd==1); K++)  {

           
      // Just for one sub-phase of size reduction
      setId(Uloc);

      if (odd==0) swapin = 0; 
    
      time.start();

      // Size reduction 
      // --------------  

    
      hsizereduce(endposition);

      if (endposition >=d)
	endposition=1;
      else 
	if ((K%2)==0) endposition*=2;

      for (j=0; j<d; j++) 
	for (i=j+1; i<d; i++)
	  R.set(i,j,0.0); 

   
      time.stop();
      tsize+=time;
  
      // Lovasz tests 
      // ------------

      for (kappa=1+odd; kappa<d; kappa+=2) {

	lovtest.mul(R.get(kappa-1,kappa-1),R.get(kappa-1,kappa-1));
	lovtest.mul(deltab,lovtest);
    
      
	nblov+=1;
    
	s.mul(R.get(kappa,kappa),R.get(kappa,kappa));
	sn.mul(R.get(kappa-1,kappa),R.get(kappa-1,kappa));
	newt.add(s,sn);

            
	if (lovtest > newt) { // Swap, down 
    
	  swapin +=1;

	  nbswaps+=1;
       

	  B.colswap(kappa-1,kappa);

	  Uloc.colswap(kappa-1,kappa);
	
	  if (transf) U.colswap(kappa-1,kappa);

	} // end swap 
     
      
      } // End column loop for Lovasz tests

      tup+=time;
     

      odd=(odd +1) % 2;
    
    
    } // Main iteration loop 


    hsizereduce(d);
  
    return 0;
  
  };
  
  
  template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
  PLattice<ZT,FT, MatrixZT, MatrixFT>::seysenreduce() { 

    householder();


    return 0;
    
  } 
  

/* -------------------------------------------------------------------------
   Size reduction 

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm

   ------------------------------------------------------------------------- */


    template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
      PLattice<ZT,FT, MatrixZT, MatrixFT>::hsizereduce(int endposition) { 

    
  FP_NR<FT> approx;
  
  approx=0.1;


  FP_NR<FT> t,tmpfp;
  Z_NR<ZT>  xz,tmpz;

  long expo,lx;

  int i,j,w=0;

  bool nonstop=1;
  //vector<bool> somedone(d);
  vector<long> somedone(d);
  
 
  householder();
 
    
  // for (j=0; j<d; j++) {
  //   col_kept[j]=0;
  //   householder_r(j);
  //   householder_v(j);
  // }

 
    
  // While loop for the norm decrease
  // --------------------------------


  //cout << "end " << endposition << endl; 
    
  while (nonstop) {

    itime.start();
    
    w++;
    
    
    for (j=1; j<d; j++) { // loop on all the columns

             
      somedone[j] = 0;

      //itime.start();
 

      // Loop through the column 
      // -----------------------

      
      //endposition=max(0,j-dec);
            
      for (i=j-1; i>=max(0,j-endposition); i--){  

         
	x.div(R.get(i,j),R.get(i,i));
	
	x.rnd(x);

	
	
	if (x.sgn() !=0) {   // Non zero combination 
	                     // --------------------
	  lx = x.get_si_exp(expo);

	  	 
	  // Cf fplll 
	  // Long case 
	  if (expo == 0) {

	    if (lx == 1) {

	      compteur +=1;
	      somedone[j] += 1;
	      
	      
	      R.subcol(j,i,i+1);

	      B.subcol(j,i,nmax);

	      //Uloc.subcol(j,i,d);
	    
	      if (transf) 
		U.subcol(j,i,d);

	    	
	    } 
	    else if (lx == -1) {

	      compteur +=1;
	      somedone[j] += 1;
	      
 
	      R.addcol(j,i,i+1);

	      B.addcol(j,i,nmax);

	      //Uloc.addcol(j,i,d);
	     
	      if (transf) 
		U.addcol(j,i,d);

	    
	    } 
	    else { 

	      compteur +=1;
	      somedone[j] += 1;
	      

	      if (fast_long_flag == 1) {
	    
		R.submulcol(j,i,x,i+1);
	       		
		B.addmulcol_si(j,i,-lx,nmax); // ICI 

		if (transf)  
		  U.addmulcol_si(j,i,-lx,d);

	      } // end fast_long
	      else {
		
		set_f(xz,x);  
		
		R.submulcol(j,i,x,i+1);	
		B.submulcol(j,i,xz,nmax);
		if (transf)  
		  U.submulcol(j,i,xz,d);
	      }

	    } // end expo == 0 and not 1 or -1
  
	  } // end expo == 0 
	  else {  // expo <> 0 

	    compteur +=1;
	    somedone[j] += 1;
	    
 
	    if (fast_long_flag == 1) {
	      
	      R.submulcol(j,i,x,i+1);

	      
	      B.addmulcol_si_2exp(j,i,-lx,expo,nmax);

	      if (transf)  
		U.addmulcol_si_2exp(j,i,-lx,expo,d);
	      
	    } // end fast_long
	    else {
	      
	      set_f(xz,x);  
	      
	      R.submulcol(j,i,x,i+1);	
	      B.submulcol(j,i,xz,nmax);
	      if (transf)  
		U.submulcol(j,i,xz,d);
	      
	    } // end no long
 
	  } // end expo <> 0 

	} // Non zero combination 

      } // Loop through the column
    
    
    } // End loop on the columns

  
    nonstop = false;

    for (j=0; j<d; j++)
      old_normB2[j]=normB2[j];

    // DBG Si somedone
     // {
     //   int nb=0;
     //   for (int k=0; k<d; k++)
     // 	if (somedone[k]>0) nb+=1;

     //   //cout << nb << " / " << d << endl;

     //   if (nb==0) cout << "****************** " << endl; 
       
     // }

    
    
    itime.stop();
    ichrono+=itime;

   
    householder();
      
    // for (j=0; j<d; j++) {
    //   if (somedone[j] >0) col_kept[j]=0;
    //   householder_r(j); 
    //   householder_v(j);
    // }

    
    
    for (j=1; j<d; j++) {
	
	t.mul(approx,old_normB2[j]);

	if (normB2[j] < t) nonstop = true;
    }
  
   
  } // end while 

  return 0;

}


/* ------------- */
/* Householder R */
/* ------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
PLattice<ZT,FT, MatrixZT, MatrixFT>::householder_r(int kappa)
{

 
  //nmaxkappa=structure[kappa]+1;
  nmaxkappa = nmax;
  
  int i,k,length; 

  // kappa == 0
  // ----------
  if (kappa==0) {

     if (col_kept[kappa])    
     
       R.setcol(kappa,Bfp.getcol(kappa),nmaxkappa);

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
    for (i=0; (i<kappa) && (i< nmaxkappa); i++) toR[i]=Rkept.get_non_normalized(i,i);

    for (i=kappa; i<nmaxkappa; i++) toR[i]=Rkept.get_non_normalized(i,kappa-1);
    
    R.setcol(kappa,&toR[0],nmaxkappa);
    
    kappamin[kappa]=kappa;
    
  } // else kappa !=0 

 return 0;
}

/* ------------- */
/* Householder V */
/* ------------- */
// nmaxkappa must be initialized e.g. through householder_r 

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
PLattice<ZT,FT, MatrixZT, MatrixFT>::householder_v(int kappa) 
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
PLattice<ZT,FT, MatrixZT, MatrixFT>::set_nblov_max(unsigned int nb) {

  nblov_max = nb;
  return nblov_max;

} 

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline unsigned int 
PLattice<ZT,FT, MatrixZT, MatrixFT>::setprec(unsigned int prec) {

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
PLattice<ZT,FT, MatrixZT, MatrixFT>::getprec() {

  return (R.get(0,0)).getprec(); 
  
}


template<class ZT,class FT,class MatrixZT, class MatrixFT> inline ZZ_mat<ZT> PLattice<ZT,FT, MatrixZT, MatrixFT>::getbase()
{
  ZZ_mat<ZT> BB(n,d);
  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) BB.Set(i,j,B(i,j)); // reprendre boucle sur les colonnes 

  return BB;
}


template<class ZT,class FT,class MatrixZT, class MatrixFT> inline MatrixZT PLattice<ZT,FT, MatrixZT, MatrixFT>::getmixedbase()
{
  return B;
}


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline ZZ_mat<ZT> PLattice<ZT,FT, MatrixZT, MatrixFT>::getU()
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

  
template<class ZT,class FT, class MatrixZT, class MatrixFT> inline  matrix<FP_NR<FT> > PLattice<ZT,FT, MatrixZT, MatrixFT>::getR()
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
PLattice<ZT,FT, MatrixZT, MatrixFT>::init(int n, int d, bool forU) {

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
  old_normB2.resize(d);

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
PLattice<ZT,FT, MatrixZT, MatrixFT>::PLattice(ZZ_mat<ZT> A, bool forU, int reduction_method) {
 
  
  n=A.getRows();
  d=A.getCols();

  init(n,d, forU); 

  int i,j;

  B.resize(n,d);  // Not in init for the mixed matrix case also 

  for (i=0; i<n; i++) 
    for (j=0; j<d; j++) 
      B(i,j)=A.Get(i,j); // Used to compute the structure below 


  if (transf) {    // Not in init for the mixed matrix case also 
    
    U.resize(d,d);
    for (i=0; i<d; i++) U(i,i)=1;     
  }

  Uloc.resize(d,d);
  
  nblov_max = 4294967295;
  
  seysen_flag=reduction_method;

  if (reduction_method == DEF_REDUCTION) fast_long_flag = 1;
  else if  (reduction_method == NO_LONG)  fast_long_flag = 0;
  
  matrix_structure(structure, B, n,d);

  nmax=structure[0];
  for (j=1; j<d; j++)
    if (nmax < structure[j]) nmax= structure[j]+1;
  
  for (j=0; j<d; j++) 
    Bfp.setcol(j,B.getcol(j),0,structure[j]+1); 


 }


// ------------------------------------------------------------------
// Assigns a basis A 
// Assumes that the lattice has already been initialized 
// 
// -------------------------------------------------------------------  

template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
PLattice<ZT,FT, MatrixZT, MatrixFT>::assign(ZZ_mat<ZT> A) {

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
PLattice<ZT,FT, MatrixZT, MatrixFT>::assign(MatrixZT A) {

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
PLattice<ZT,FT, MatrixZT, MatrixFT>::shift_assign(ZZ_mat<ZT> A, vector<int> shift) {

  nblov = 0;

  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) 
      B(i,j).mul_2si(A(i,j),shift[i]);
  
  
  // ICI TMP 
  /*
  ZZ_mat<ZT> T;
  T.resize(n,d);

  for (int i=0; i<n; i++) 
    for (int j=0; j<d; j++) 
      T(i,j)=B(i,j);
      
  int lnorm=0;
  Z_NR<ZT> norm;

  
  for(int j=0; j<d; j++)  {
    
    norm.mul(T(0,j),T(0,j));
    for (int i=1; i<n; i++)
      norm.addmul(T(i,j),T(i,j));
    
    lnorm = max(lnorm,size_in_bits(norm)/2);
	    
    }
  
  //cout << "lnorm " << " : " << lnorm << endl; 
  
  for (int i=0; i<n; i++)
    for (int j=0; j<d; j++) 
      T(i,j).mul_2si(T(i,j),-lnorm+max(d,40)+8);

  //print2maple(T,n,d);
  //print2maple(B,n,d);

  for (int i=0; i<n; i++)
    for (int j=0; j<d; j++) 
    B(i,j)=T(i,j);
  
  //print2maple(T,n,d);
  */


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
PLattice<ZT,FT, MatrixZT, MatrixFT>::put(ZZ_mat<ZT> A, long upperdim, long t, long tau) {

  // Cf structure, colkept, kappamin 

  trunc<ZT, MatrixZT>(B, A, upperdim, d, n, t, tau);

  if (transf) {

    U.resize(d,d);
    for (int i=0; i<d; i++) U(i,i)=1; 

  }
}


template<class ZT,class FT, class MatrixZT, class MatrixFT>  void 
PLattice<ZT,FT, MatrixZT, MatrixFT>::mixed_put(MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> > A, long t, long tau) {

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
PLattice<ZT,FT, MatrixZT, MatrixFT>::shift(ZZ_mat<ZT> A, long m, long lsigma) {
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
PLattice<ZT,FT, MatrixZT, MatrixFT>::shiftRT(long lsigma) {
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

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline void PLattice<ZT,FT, MatrixZT, MatrixFT>::isreduced(double deltain) {

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
  cout << endl;
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
PLattice<ZT,FT, MatrixZT, MatrixFT>::lcond(int tproper, int flagprec,  int flagheur) {


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

    PLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(getbase());

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

      PLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(getbase());
      
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
    
    PLattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > L(getbase());
    
    return(L.lcond(ANY, DEFAULT_PREC));
    
    mpfr_set_default_prec(oldprec);


  } // and unknown 

  return 0; //cc
}




/* --------------------------------------------- */
/* Householder complet */
/* --------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
PLattice<ZT,FT, MatrixZT, MatrixFT>::householder()
{

  int i,k,kappa;
  FP_NR<FT> nrtmp,s,w; 
  
  nrtmp=0;
  s=0;

  

    for (kappa=0; kappa<d; kappa++) {

      R.setcol(kappa,B.getcol(kappa),0,nmax);

      fp_norm_sq(normB2[kappa], R.getcol(kappa), nmax);
       
      for (k=0; k<kappa; k++) {
	scalarprod(nrtmp, V.getcol(k,k), R.getcol(kappa,k), nmax-k);
	R.fmasub(kappa,k,R.getcol(kappa,k), V.getcol(k,k), nrtmp, nmax-k); 
      }


      w=R.get(kappa,kappa);

      if (w >=0) {
	fp_norm(s,R.getcol(kappa,kappa),nmax-kappa); 
	nrtmp.neg(s);
	R.set(kappa,kappa,nrtmp);    
      }
      else {
	fp_norm(nrtmp,R.getcol(kappa,kappa),nmax-kappa); // de la colonne
	R.set(kappa,kappa,nrtmp);
	s.neg(nrtmp);  
      }

      w.add(w,s);
      s.mul(s,w);
      s.sqrt(s);

      V.div(kappa,kappa+1, R.getcol(kappa,kappa+1), s, nmax-kappa-1);

      nrtmp.div(w,s);
      V.set(kappa,kappa,nrtmp); 

      for(i=kappa+1; i<d; i++)  R.set(i,kappa,0.0); 
       
    }  // sur kappa 
   
    return 0; 
}




/* -------------------------------------------------------------------------
   Size reduction 

   Assumes that Householder is available until index kappa-1

   Returns -1 if no convergence in the while loop:  no decrease of the norm

   ------------------------------------------------------------------------- */



// The swap is done already 
template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
PLattice<ZT,FT, MatrixZT, MatrixFT>::qrupdate(int iend) { 

  int kappa,j;

  FP_NR<FT> t1,t2,x,q11,q22,q12,q21; 

  // QR update 
  // ---------

  for (kappa=1; kappa <d ; kappa++) { 

    if ((R.get(kappa,kappa-1)).sgn() !=0) { // Not diagonal 
    
	
      t1.mul(R.get(kappa-1,kappa-1),R.get(kappa-1,kappa-1));
      t2.mul(R.get(kappa,kappa-1),R.get(kappa,kappa-1));
      x.add(t1,t2);
      x.sqrt(x);

      q11.div(R.get(kappa-1,kappa-1),x);
      q22=q11;
      q12.div(R.get(kappa,kappa-1),x);
      q21.neg(q12);

      // Col kappa-1
      R.set(kappa-1,kappa-1,x);
      R.set(kappa,kappa-1,0.0);

      // Other columns, could be restricted to only some columns
      // -------------------------------------------------------


      for (j=kappa; j<d; j++) {

	t1.mul(q11,R.get(kappa-1,j));
	t2.mul(q12,R.get(kappa,j));
	t1.add(t1,t2);

	
	t1.mul(q21,R.get(kappa-1,j));

	R.set(kappa-1,j,t1); // after it is used for next 

	t2.mul(q22,R.get(kappa,j));
	t1.add(t1,t2);

	R.set(kappa,j,t1);
      }       
    
    } // end if non zero for diag 
  } // end loop on the columns 

  // Size reduction 
  // --------------

  FP_NR<FT> approx;
  approx=0.1;

  FP_NR<FT> t,tmpfp;
  Z_NR<ZT>  xz,tmpz;
  
  long expo,lx;
  
  int i;
  
  int nnmax; // De la structure triangulaire 
  
  for (kappa=1 ; kappa<d; kappa++) { 
    
    nmaxkappa=structure[kappa]+1;
  

    // Loop through the column 
    // -----------------------

    int limit;
   
    if (iend ==1) limit =-1;
    else limit = kappa-2;

    for (i=kappa-1; i>limit; i--){  
         
      x.div(R.get(i,kappa),R.get(i,i)); 
      x.rnd(x);

      if (x.sgn() !=0) {   // Non zero combination 
                           // --------------------
	lx = x.get_si_exp(expo);

	nnmax=structure[i]+1;
	
	// Cf fplll 
	// Long case 
	if (expo == 0) {

	  if (lx == 1) {

	    R.subcol(kappa,i,i+1);	    
	    Bfp.subcol(kappa,i,n);
      	    if (transf) 
	      U.subcol(kappa,i,min(d,nnmax));
	    
	  } 
	  else if (lx == -1) {
 
	    R.addcol(kappa,i,i+1);
	    Bfp.addcol(kappa,i,n);
	    if (transf) 
	      U.addcol(kappa,i,min(d,nnmax));
	    

	  } 
	  else { 
 
	    R.submulcol(kappa,i,x,i+1);
	    Bfp.submulcol(kappa,i,x,n);
	    if (transf) 
	      U.addmulcol_si(kappa,i,-lx,min(d,nnmax));
	    
	  } 
  
	} // end expo == 0 
	else {  // expo <> 0 

	  set_f(xz,x);
	  R.submulcol(kappa,i,x,i+1);
	  //B.submulcol(kappa,i,xz,nnmax);
	  Bfp.submulcol(kappa,i,x,n);
	  if (transf)  
	    //U.submulcol(kappa,i,xz,min(d,nnmax));
	    U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nnmax));
	  
	} // end expo <> 0 
      } // Non zero combination 

    } // Loop through the     
  } // end loop over the columns     
 

  //==============
  return 0;  
} 


} // end namespace hplll

#endif 

