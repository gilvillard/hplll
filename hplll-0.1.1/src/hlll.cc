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


#ifndef HPLLL_HLLL_CC
#define HPLLL_HLLL_CC


// MPFR setprec : isreduced and cond are re-initializing the floating point matrices 

namespace hplll { 

template<class ZT,class FT, class MatrixZT, class MatrixFT>  int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::hlll(double delta, bool verbose) { 

  int kappa=1,i;
  int prevkappa=-1; // For the looping test betwenn tow indices 
  vector<FP_NR<FT> >  prevR(d);


  FP_NR<FT> newt; //testaccu;

  FP_NR<FT> deltab,lovtest;
  deltab=delta;   // TO SEE 

  FP_NR<FT> tmpswap;

  FP_NR<FT> s,sn; // newt test 

  unsigned int start=0,starttot;
  starttot=utimesec();

  int flag_reduce=0; // No convergence in reduce 


  for (i=0; i<d; i++) {col_kept[i]=0; descendu[i]=0;}

  // ICI 
  

  while ((kappa < d) && (nblov < nblov_max)) 
    {

      if (((nblov%800000)==0) && (nblov > 0))   cout << nblov << " tests" << endl; 
      
      if (kappa == 1) { 
	for (i=0; i<d; i++) kappamin[i]=min(kappamin[i],0);
	householder_r(0);  
	householder_v(0);
	
      } 
      else for (i=kappa+1; i<d; i++) kappamin[i]=min(kappamin[i],kappa); // lowest visit after kappa 

      
      if (chrono)  start=utime();
      
     
      // ICI modif descendu 
      if (descendu[kappa]>=1) {
	
	if (seysen_flag == 0)
	  flag_reduce=hsizereduce(kappa);   
	else 
	  flag_reduce=seysenreduce(kappa); 
	
      }
      else { 
	
	if (seysen_flag == 0)
	  flag_reduce=hsizereduce(kappa);   
	else 
	  flag_reduce=seysenreduce(kappa); 
	
      }
      
      if (chrono) tps_reduce+=utime()-start;
      
      if (flag_reduce==-1) return(-1);
      
      if (chrono) start=utime();

      //newt=normB2[kappa];  // The newt test was like that in old non exp 
      //for (i=0; i<=kappa-2; i++) newt.submul(R.get(i,kappa),R.get(i,kappa));

      lovtest.mul(R.get(kappa-1,kappa-1),R.get(kappa-1,kappa-1));
      lovtest.mul(deltab,lovtest);

      if (chrono) tps_prepare+=utime()-start;

      if (chrono) start=utime();
   
      nblov+=1;
    
      fp_norm(s,R.getcol(kappa,kappa),structure[kappa]+1-kappa);
      s.mul(s,s);
      sn.mul(R.get(kappa-1,kappa),R.get(kappa-1,kappa));
      newt.add(s,sn);
 

      // ****************
      //   UP    UP    UP 
      // ****************
    
      if (lovtest <= newt) {
	
	if (verbose) {
	  if (((kappa) < d) && (col_kept[kappa+1] == 0)) cout << "Discovering vector " 
							      << kappa +1 << "/" << d 
							      << " cputime=" << utimesec() -starttot << "sec" << endl; 
	}

	householder_v(kappa);   // The first part of the orthogonalization is available 
	if (chrono) tps_householder+=utime()-start;

	

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
	if (lsize > 0) L.colswap(kappa-1,kappa);

	Bfp.colswap(kappa-1,kappa);

	structure[kappa-1]=structure[kappa];

	VR.colswap(kappa-1,kappa);
 
	tmpswap=normB2[kappa];  normB2[kappa]=normB2[kappa-1]; normB2[kappa-1]=tmpswap;

	prevkappa=kappa; 
	kappa=max(kappa-1,1);  

      }

    if (chrono) tps_swap+=utime()-start;    
   

    } // End main LLL loop 
  
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


  FP_NR<FT> x,t,tmpfp;
  Z_NR<ZT>  xz,tmpz;

  long expo,lx;

  vector<FP_NR<FT> > vectx(kappa);  
  
  vector<FP_NR<FT> > tmpcolR(kappa);  

  int i,k,w=0;

  bool nonstop=1;
  bool somedone=0;

  unsigned int start=0;
 
  int restdim=0; // Remaining dimension after the current block 

  int nmax; // De la structure triangulaire 

  //int whilemax=10000;  // Convergence problem if too low precision  

  // To see / prec problem 
  //col_kept[kappa]=0;

  if (chrono)  start=utime();
  householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 

  if (chrono) tps_householder+=utime()-start;

  // DBG to SEE 
  //col_kept[kappa]=0;

  int bdim,ld,tdig,indexdec;

 
  while (nonstop) {  // LOOP COLUMN CONVERGENCE

    
    w++;

    somedone = 0;

    if (chrono) start=utime();
    
    // ----------  NOUVELLE BOUCLE, SUR LES BLOCS 
    // Kappa est la dimension de ce qu'il y a avant 
     
    //DBG          cout << "********** kappa  : " << kappa << endl; 

    ld=1; indexdec=0; // Décalage d'indice
    while (ld <=kappa) {

      tdig=(kappa/ld)%2;
      if (tdig==0) bdim =0; else bdim = ld;

      // -----------------------------------------------------
      // Boucle sur la partie de la colonne correspond au bloc 
      // -----------------------------------------------------

      // vectxz rouding of vectx
      //  column -  (prev col) * vectxz (xz rÃ©-utilisÃ©) 
      // On peut travailler sur place en remontant dans la colonne kappa de R 

      // On calcule vectx et on arrondit au fur et Ã  mesure

      restdim=kappa-indexdec-bdim;

      //DBG 
      //householder_r(kappa);

      for (i=kappa-1-indexdec; i>=restdim; i--) 
	tmpcolR[i]=R.get(i,kappa);

      for (i=kappa-1-indexdec; i>=restdim; i--){
	 
	vectx[i].div(tmpcolR[i],R.get(i,i));
	for (k=restdim; k<i; k++) tmpcolR[k].submul(R.get(k,i),vectx[i]);

	vectx[i].rnd(vectx[i]);

      } // end calcul de la transfo 

      // Et on applique la transformation  
      // --------------------------------
      for (i=kappa-1-indexdec; i>= restdim; i--){
    
	x=vectx[i]; 

	if (x.sgn() !=0) { 

	  lx = x.get_si_exp(expo);

	  nmax=structure[i]+1;
	  
	  // Cf fplll 
	  // Long case 
	  if (expo == 0) {
	    
	    if (lx == 1) {
	      
	      somedone = 1;
	      
	      R.subcol(kappa,i,restdim);
	      
	      B.subcol(kappa,i,nmax);
	      
	      if (transf) 
		U.subcol(kappa,i,min(d,nmax));

	      if (lsize > 0) 
		L.subcol(kappa,i,lsize);
	      
	    } 
	    else if (lx == -1) {

	      somedone = 1;
 
	      R.addcol(kappa,i,restdim);
	      
	      B.addcol(kappa,i,nmax);
	      
	      if (transf) 
		U.addcol(kappa,i,min(d,nmax));

	      if (lsize > 0) 
		L.addcol(kappa,i,lsize);

	    } 
	    else { 
 
	      somedone = 1;
	      
	      R.submulcol(kappa,i,x,restdim);
	      
	      B.addmulcol_si(kappa,i,-lx,nmax);
	      
	      if (transf) 
		U.addmulcol_si(kappa,i,-lx,min(d,nmax));

	      if (lsize > 0) 
		L.addmulcol_si(kappa,i,-lx,lsize);
	    } 
	    
	  } // end expo == 0 
	  else {  // expo <> 0 

	    somedone = 1;

	    set_f(xz,x);
	    
	    R.submulcol(kappa,i,x,restdim);

	    if (chrono) start=utime();
	     
	    B.addmulcol_si_2exp(kappa,i,-lx,expo,nmax);

	    if (transf)  
	      U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));

	    if (lsize >0)  
	      L.addmulcol_si_2exp(kappa,i,-lx,expo,lsize);

	    if (chrono) tps_redB+=utime()-start;	 
	  } // end expo <> 0 
	} // Non zero combination 

      } // end application de la transformation 

      indexdec+=bdim;     
      ld=ld*2;

    } // End loop on log blocks 

    if (chrono)  tps_redB+=utime()-start;

    

    if (somedone) {
     
      col_kept[kappa]=0;

      t.mul(approx,normB2[kappa]);

      if (chrono) start=utime();
      householder_r(kappa);  
      if (chrono) tps_householder+=utime()-start;

      nonstop = (normB2[kappa] < t);  // ne baisse quasiment plus ? 

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
Lattice<ZT,FT, MatrixZT, MatrixFT>::hsizereduce(int kappa) { 

  nmaxkappa=structure[kappa]+1;

  FP_NR<FT> approx;
  
  approx=0.1;


  FP_NR<FT> t,tmpfp;
  Z_NR<ZT>  xz,tmpz;

  long expo,lx;

  int i,w=0;

  bool nonstop=1;
  bool somedone=0;

  unsigned int start=0;

  int nmax; // De la structure triangulaire 
  
  if (chrono) start=utime();
  
  householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 
  
  if (chrono) tps_householder+=utime()-start;


  // While loop for the norm decrease
  // --------------------------------
  

  while (nonstop) {

    w++;

    somedone = 0;

    if (chrono) start=utime();

    // Loop through the column 
    // -----------------------

    for (i=kappa-1; i>-1; i--){  

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

	    somedone = 1;

	    R.subcol(kappa,i,i+1);
	    
	    B.subcol(kappa,i,nmax);
		
	    if (transf) 
	      U.subcol(kappa,i,min(d,nmax));

	    if (lsize > 0) 
	      L.subcol(kappa,i,lsize);
	
	  } 
	  else if (lx == -1) {

	    somedone = 1;
 
	    R.addcol(kappa,i,i+1);
	    
	    B.addcol(kappa,i,nmax);
		
	    if (transf) 
	      U.addcol(kappa,i,min(d,nmax));

	    if (lsize > 0) 
	      L.addcol(kappa,i,lsize);

	  } 
	  else { 
 
	    somedone = 1;

	    R.submulcol(kappa,i,x,i+1);
	   
	    B.addmulcol_si(kappa,i,-lx,nmax);


	    if (transf) 
	      U.addmulcol_si(kappa,i,-lx,min(d,nmax));

	    if (lsize >0) 
	      L.addmulcol_si(kappa,i,-lx,lsize);

	  } 
  
	} // end expo == 0 
	else {  // expo <> 0 

	  somedone = 1;

	  set_f(xz,x);
	  
	  R.submulcol(kappa,i,x,i+1);
      

	  B.submulcol(kappa,i,xz,nmax);    
	  //B.addmulcol_si_2exp(kappa,i,-lx,expo,nmax);
	 
	  if (transf)  
	    U.submulcol(kappa,i,xz,min(d,nmax));
	  //U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));

	  if (lsize > 0)  
	    L.addmulcol_si_2exp(kappa,i,-lx,expo,lsize);

	} // end expo <> 0 

      } // Non zero combination 

    } // Loop through the column
    
    if (chrono)  tps_redB+=utime()-start;

  

    if (somedone) {
 
      col_kept[kappa]=0;

      t.mul(approx,normB2[kappa]);
      householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 
      nonstop = (normB2[kappa] < t);  // ne baisse quasiment plus ? 
      
    }
    else {
      nonstop=0;
      
    }
  } // end while 
   
 

  return somedone;

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
	  
	  for (k=1; k<kappa; k++) {
	    length--;
	    scalarprod(VR(k,kappa), V.getcol(k,k), Rkept.getcol(k-1,k), length);
	    Rkept.fmasub(k,k, Rkept.getcol(k-1,k), V.getcol(k,k), VR(k,kappa), length);  //  k-1 to k 
	  }
	}
	else {
	  
	  k=0;
	  
	  Rkept.fmasub(k,k, Bfp.getcol(kappa,k), V.getcol(k,k), VR(k,kappa), length);  // k=0
	  
	  for (k=1; k<kappamin[kappa]; k++)  {
	    length--;
	    Rkept.fmasub(k,k, Rkept.getcol(k-1,k), V.getcol(k,k), VR(k,kappa), length);  // k-1 to k
	  } 
	  
	  for (k=kappamin[kappa]; k<kappa; k++) {
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
      
      for (k=1; k<kappa; k++) {
	
	length--;
	scalarprod(VR(k,kappa), V.getcol(k,k), Rkept.getcol(k-1,k), length);
	
	Rkept.fmasub(k,k, Rkept.getcol(k-1,k), V.getcol(k,k), VR(k,kappa), length);  // de k-1 à k 
      }
      

      length = nmaxkappa; 

    } // endelse recomputation 
    // -----------------------------------------------------------
    
     
    // Dummy in the standard case, made special for the MatrixPE case 
   
    for (i=0; i<kappa; i++) toR[i]=Rkept.get_non_normalized(i,i);
    
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
Lattice<ZT,FT, MatrixZT, MatrixFT>::householder_v(int kappa) 
{
  int i;
  FP_NR<FT> s,norm,w,tmpdpe; 

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

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline ZZ_mat<ZT> Lattice<ZT,FT, MatrixZT, MatrixFT>::getL()
{
 
  if (lsize > 0) { 
    ZZ_mat<ZT> LL(lsize,d); 
    for (int i=0; i<lsize; i++) 
      for (int j=0; j<d; j++) LL.Set(i,j,L.get(i,j)); // reprendre boucle sur les colonnes 
    return LL;
  }
  else {
    cout << "*** Error, HLLL, the Lehmer companion matrix has not been computed" << endl;
    ZZ_mat<ZT> LL(0,0);
    return LL;
  }
}


template<class ZT,class FT, class MatrixZT, class MatrixFT> inline  matrix<FP_NR<FT> > Lattice<ZT,FT, MatrixZT, MatrixFT>::getR()
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
Lattice<ZT,FT, MatrixZT, MatrixFT>::Lattice(ZZ_mat<ZT> A, bool forU, int reduction_method, int lehmer_size) {

  
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

  nblov_max = 4294967295;
  lsize = 0; // For other cases 

  if (lehmer_size >0) {

    lsize = lehmer_size;
    L.resize(lsize,d);

    for (i=0; i<lsize; i++) 
      for (j=0; j<d; j++) 
	L(i,j)=B(i,j);

  }

  seysen_flag=reduction_method;

  matrix_structure(structure, B, n,d);


 }


// With mixed matrix from mixed matrix 
// -----------------------------------

template<class ZT,class FT, class MatrixZT, class MatrixFT>
Lattice<ZT,FT, MatrixZT, MatrixFT>::Lattice(MatrixRZ<matrix, FP_NR<mpfr_t>, Z_NR<ZT> > A, bool forU, int reduction_method) {

    n=A.getRowsRT()+A.getRowsZT();
    d=A.getCols();

    init(n,d, forU); 
    
    int i,j;

    B.resize(A.getRowsRT(), A.getRowsZT(), d);

    for (i=0; i<A.getRowsRT(); i++) 
      for (j=0; j<d; j++) 
	B.setRT(i,j,A.getRT(i,j));

    for (i=0; i<A.getRowsZT(); i++) 
      for (j=0; j<d; j++) 
	B.setZT(i,j,A.getZT(i,j));
    
    if (transf) {    // Not in init for the mixed matrix case also 
      Z_NR<ZT> one;
      one=1;
      U.resize(0,d,d);
      for (i=0; i<d; i++) U.set(i,i,one); 
    }

    nblov_max = 4294967295;
    lsize = 0; // For other cases 

    seysen_flag=reduction_method;
  
    matrix_structure(structure, B, A.getRowsRT(), A.getRowsZT(), d);
}


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
Lattice<ZT,FT, MatrixZT, MatrixFT>::shift_assign(ZZ_mat<ZT> A, vector<int> shift) {

  nblov = 0;

  if (lsize > 0) {

    for (int i=0; i<lsize; i++) 
      for (int j=0; j<d; j++) 
	B(i,j).mul_2si(L(i,j),shift[i]);


  }
  else {
    for (int i=0; i<n; i++) 
      for (int j=0; j<d; j++) 
	B(i,j).mul_2si(A(i,j),shift[i]);
  }
  
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

  if (nblov==d-1) cout << "Seems to be reduced (up to size-reduction)" << endl << endl;

  if (nblov != d-1) {
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
// Condition number || |R| |R^-1| || _ F
// Change the precision, hence initialize R (forget previous value)  

// To be used with mpfr if the precision change must have some effect 
//*************************************************************************** 
// Matrix interface () deprecated, will note work with -exp 
//

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline FP_NR<FT> Lattice<ZT,FT, MatrixZT, MatrixFT>::cond() {


  FP_NR<FT>  cc,c2;
 

  int i,j,k;

  unsigned int oldprec,newprec;

  oldprec=getprec();

  // ICI 
  //if (d<=20) setprec(53);  
  //else setprec(2*d);  // ******* TO TUNE  

  newprec=getprec();

  //if (newprec == oldprec) cout << "Warning: in function cond the change of precision has no effect" << endl; 

  // Calcul direct sur R par householder 
  // -----------------------------------

  householder();

  // Abs R --> ar 
  // Do not change for the template class MatrixZT, class MatrixFT that is different e.g. double and mpfr 
  MatrixFT aR;
  aR.resize(d,d);

  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      aR(i,j).abs(R(i,j));

  // Magnitude of the diagonal 
  // -------------------------

  FP_NR<FT>  mindiag,maxdiag;
  mindiag=aR(0,0);
  maxdiag=aR(0,0);

  for (i=1; i<d; i++) {
    if (mindiag > aR(i,i)) mindiag=aR(i,i);
    if (maxdiag < aR(i,i)) maxdiag=aR(i,i);
  }
  
  //cout << endl << "Diagonal max: " << maxdiag << "  min: " << mindiag << endl;
  maxdiag.div(maxdiag,mindiag);
  //cout << "Largest diagonal ratio: " << maxdiag << endl;

  // Inversion 
  // ---------


  MatrixFT iR;
  iR.resize(d,d);

  for (i=0; i<d; i++) iR(i,i)=1.0; 

  //  for (i=1; i<=n; i++) SET_D(II[i][i], 1.0);

  // Higham p263 Method 1

  FP_NR<FT> t,one;
  one=1.0;

  for (j=0; j<d; j++)
    {

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


  cc=0.0; 

  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      cc.addmul(prod(i,j),prod(i,j));

  cc.sqrt(cc); 

  cout << endl;
  cout << "Conditioning is about ";
  cc.print();
  cout  << endl; 

  cc.log(cc);
  c2=2.0;
  c2.log(c2);
  cc.div(cc,c2);
  cout << "Digits - Log conditioning is about ";
  cc.print();

  Z_NR<long> tt;
  tt.set_f(cc);
  cout << " = " << tt << endl; 
  cout  << endl; 


  setprec(oldprec);

  return cc;
}



//***************************************************************************
// Version simplifiée de l'énergie
// Change la prec donc écrase R 

// Avec mpfr si le chgt de précision doit avoir un effet 
// Matrix interface () deprecated, will note work with -exp 
//

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline FP_NR<FT> Lattice<ZT,FT, MatrixZT, MatrixFT>::energy() {

  
  unsigned int oldprec,newprec;

  oldprec=getprec();

  if (d<=20) setprec(53);  
  else setprec(2*d);  // ******* TO TUNE  

  newprec=getprec();

  if (newprec == oldprec) cout << "Warning: in function energy, the change of precision has no effect" << endl; 


  // Calcul direct sur R par householder 
  // -----------------------------------

  int i,k;

  
  householder();

  FP_NR<FT>  cc;

  FP_NR<FT> tmp,nf;

  cc=0.0;
  
  for (i=0; i<d; i++) {
    tmp.abs(R(i,i)); 
    //cout << R(i,i) << endl; 
    tmp.log(tmp);
    nf=0.0;
    for (k=0; k<d-i; k++)  // multiplication par i ! 
      nf.add(nf,tmp); 
    cc.add(cc,nf);
    }

  cout << endl;
  cout << "Energy is about ";
  cc.print();
  cout  << endl; 


  setprec(oldprec);

  return cc;
}




/* --------------------------------------------- */
/* Householder complet */
/* --------------------------------------------- */

template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::householder()
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



//***************************************************************************
// Stats for the nullspace tests 

// Avec mpfr si le chgt de precision doit avoir un effet 
// Time in input for file printing 
/*
template<class ZT,class FT, class MatrixZT, class MatrixFT> inline FP_NR<FT> Lattice<ZT,FT, MatrixZT, MatrixFT>::tnull(int time) {

  FP_NR<FT>  cc,c2;
 

  int i,j,k;

  int oldprec;

  if (d<=20) oldprec=setprec(53);  
  else oldprec=setprec(20*d);  // ******* A VOIR 


  // Calcul direct sur R par householder 
  // -----------------------------------

  householder();

  // Abs R --> ar 
  MatrixFT aR;
  aR.resize(d,d);

  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      aR(i,j).abs(R(i,j));

  // Magnitude of the diagonal 
  // -------------------------

  FP_NR<FT>  mindiag,maxdiag;
  mindiag=aR(0,0);
  maxdiag=aR(0,0);

  for (i=1; i<d; i++) {
    if (mindiag > aR(i,i)) mindiag=aR(i,i);
    if (maxdiag < aR(i,i)) maxdiag=aR(i,i);
  }
  
  //  cout << endl << "Diagonal max: " << maxdiag << "  min: " << mindiag << endl;
  maxdiag.div(maxdiag,mindiag);

  c2=2.0;
  c2.log(c2);
  maxdiag.log(maxdiag); 
  maxdiag.div(maxdiag,c2);
  Z_NR<ZT> rdiag;
  set_f(rdiag,maxdiag);

  //cout << "Largest diagonal ratio: " << rdiag << endl;

  // Inversion 
  // ---------


  MatrixFT iR;
  iR.resize(d,d);

  for (i=0; i<d; i++) iR(i,i)=1.0; 

  //  for (i=1; i<=n; i++) SET_D(II[i][i], 1.0);

  // Higham p263 Method 1

  FP_NR<FT> t,one;
  one=1.0;

  for (j=0; j<d; j++)
    {

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


  cc=0.0; 

  for (i=0; i<d; i++)
    for (j=0; j<d; j++) 
      cc.addmul(prod(i,j),prod(i,j));

  cc.sqrt(cc); 

  // Cond 
  // ----

  cout << endl;
  cc.log(cc);
  c2=2.0;
  c2.log(c2);
  cc.div(cc,c2);
  cout << "**** Output file stats *****" << endl << endl;
  //cout << "Digits - Log conditioning is about " << cc << endl << endl; 
  Z_NR<ZT> rcond;
  set_f(rcond,cc);
  //cout << "Rounded log cond " << rcond << endl;

  // Maxbitsize 
  // ----------

  int maxbit=0;

  int nn=B.getRows();
  int dd=B.getCols();

  for (int ii=0; ii<nn ; ii++) 
    for (int jj=0; jj<dd; jj++)
      maxbit=max(maxbit, (int) mpz_sizeinbase(B(ii,jj).getData(),2)); 

  //cout << "Max bit size " << maxbit << endl;

  // b_1
  // ---

  int maxbit1=0;
  for (int ii=0; ii<nn ; ii++) 
    maxbit1=max(maxbit1, (int) mpz_sizeinbase(B(ii,0).getData(),2)); 

  //cout << "b_1 bit size " << maxbit1 << endl;


  // With Vol 

  int kk;
  kk=B.getCols();
  
  //cout << "k  " << kk << endl;

  // vol 
  FP_NR<mpfr_t> s;
  s=R(0,0);
  for (k=1; k<d; k++) 
    s.mul(s,R(k,k));
  s.abs(s);
  //cout << " *********** " << "Vol(L) is about " << s << endl; 

  cc.log(s);
  c2=2.0;
  c2.log(c2);
  cc.div(cc,c2);

  //cout << " *********** " << "log Vol(L) is about " << cc << endl; 

  Z_NR<ZT> tmpk;
  tmpk=kk;
  FP_NR<mpfr_t> kkf;
  set_z(kkf,tmpk);
  FP_NR<mpfr_t> expected;
  expected.div(cc,kkf);
  //cout << " *********** " << "log Vol(L)1/k is about " << expected << endl; 

  kkf.sqrt(kkf);
  cc.log(kkf);
  c2=2.0;
  c2.log(c2);
  cc.div(cc,c2);

  expected.add(expected,cc);

  //cout << " *********** " << "expected  " << expected << endl; 

  Z_NR<ZT> rlambda;
  set_f(rlambda,expected);
  //cout << "Rounded expected bits  " << rlambda << endl;

  Z_NR<ZT> delta;
  delta=maxbit;
  delta.sub(delta,rlambda);

  // File editing: append to out.txt 
  // -------------------------------

  ofstream outfile; 

  outfile.open ("out.txt",ios_base::app);
  outfile << "---- " << B.getRows()-B.getCols() << " x " << B.getRows() << " / " <<
    B.getRows()<< " x " << B.getCols() <<endl; 
  outfile << "*l_max: " << maxbit << "    *l_b1: " << maxbit1 << "    *l_cond: " << rcond << "    *l_diag: " << rdiag << "    *l_Gauss: " << rlambda  
	  << "    *mdelta bits: " << delta << "   Time: " << time << " sec." << endl << endl;    
  outfile.close();

  cout << "---- " << B.getRows()-B.getCols() << " x " << B.getRows() << " / " <<
    B.getRows()<< " x " << B.getCols() <<endl; 
  cout << "*l_max: " << maxbit << "    *l_b1: " << maxbit1 << "    *l_cond: " << rcond << "    *l_diag: " << rdiag << "    *l_Gauss: " << rlambda  
       << "    *mdelta bits: " << delta <<  "   Time: " << time << " sec." << endl << endl;   

  // Precision 
  setprec(oldprec);

  return cc;
}
*/
} // end namespace hplll

#endif 

