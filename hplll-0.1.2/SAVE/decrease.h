





template<class ZT,class FT, class MatrixZT, class MatrixFT> inline int 
Lattice<ZT,FT, MatrixZT, MatrixFT>::decrease(int kappa) { 

  FP_NR<FT> approx;
  
  approx=0.01;


  FP_NR<FT> t,tmpfp;
  Z_NR<ZT>  xz,tmpz;

  long expo,expoAdd,lx;

  int i,w=0;

  bool nonstop=1;
  bool somedone=0;

  unsigned int start=0;

  int nmax; // De la structure triangulaire 

  if (chrono) start=utime();
  
  householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 

  if (chrono) tps_householder+=utime()-start;

  // ICI 
  //if (descendu[kappa] >=1) nonstop=0;

  // While loop for the norm decrease
  // --------------------------------
  // ICI 

  while (nonstop) {
 
    w++;

    somedone = 0;

    if (chrono) start=utime();

    // Loop through the column 
    // -----------------------
    
    for (i=kappa-1; i>kappa-2; i--){  
                               

      x.div(R.get(i,kappa),R.get(i,i)); 
      //x.div(R.get_non_normalized(i,kappa),R.get_non_normalized(i,i)); 
      x.rnd(x);

      
      if (x.sgn() !=0) {   // Non zero combination 
                           // --------------------
	lx = x.get_si_exp_we(expo, expoAdd);

	if (descendu[kappa] >=1) cout << "++++ " <<  lx << endl; 

	nmax=structure[i]+1;

	// Cf fplll 
	// Long case 
	if (expo == 0) {

	  if (lx == 1) {

	    somedone = 1;
	    
	    R.subcol(kappa,i,i);
	    
	    B.subcol(kappa,i,nmax);
	    
	    if (transf) 
	      U.subcol(kappa,i,min(d,nmax));
	    
	  } 
	  else if (lx == -1) {

	    somedone = 1;
	    
	    R.addcol(kappa,i,i);
	    
	    B.addcol(kappa,i,nmax);
	    
	    if (transf) 
	      U.addcol(kappa,i,min(d,nmax));
	    
	  } 
	  else { 
 
	    somedone = 1;

	    R.submulcol(kappa,i,x,i);
	   
	    B.addmulcol_si(kappa,i,-lx,nmax);


	    if (transf) 
	      U.addmulcol_si(kappa,i,-lx,min(d,nmax));
	  } 
	} // end expo == 0 
	else {  // expo <> 0 

	 

	  somedone = 1;

	  set_f(xz,x);
	  
	  R.submulcol(kappa,i,x,i);
      

	  //B.submulcol(kappa,i,xz,nmax);    
	  B.addmulcol_si_2exp(kappa,i,-lx,expo,nmax);
	 
	  if (transf)  
	    //U.submulcol(kappa,i,xz,min(d,nmax));
	    U.addmulcol_si_2exp(kappa,i,-lx,expo,min(d,nmax));

	} // end expo <> 0 

      } // Non zero combination 

    } // Loop through the column
    
    if (chrono)  tps_redB+=utime()-start;

    if (somedone) {

      col_kept[kappa]=0;

      t.mul(approx,normB2[kappa]);

      if (chrono) start=utime();
      
      householder_r(kappa); // pas tout householder necessaire en fait cf ci-dessous 
      
      if (chrono) tps_householder+=utime()-start;

      nonstop = (normB2[kappa] < t);  // ne baisse quasiment plus ? 
      //if (w > 1) cout << "----------" << endl; 
    }
    else 
      nonstop=0;

  } // end while 


  x.div(R.get(kappa-1,kappa),R.get(kappa-1,kappa-1)); 
  x.abs(x);
  if (x.cmp(0.5) > 0)  {

    cout << "++++ " << x << "   "  << kappa+1 << endl; 
    print2maple(getbase(),n,d);
    print2maple(getR(),d,d);

  }

  compteur = max(compteur,w);
   
  return somedone;

}

