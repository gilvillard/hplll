
#include <hplll.h>


// *************************************************************
// requires output matrix C and vector delta to be allocated
// delta : ith denominators (minors), starting from 1
//
// todo: to d+1? 
// *************************************************************

int ff(ZZ_mat<mpz_t>& C, vector< Z_NR<mpz_t> >& delta, ZZ_mat<mpz_t> A) {

  int i,j,k;
  
  int d=A.getCols();

  matprod(C,transpose(A),A);
  
  delta[0]=1;

  Z_NR<mpz_t> tz;
  
  for (i=0; i<d-1; i++) { // row that eliminates 

    for (k=i+1; k<d; k++) {

      // Row k with row i
      for (j=i+1; j<d; j++) {
	C(k,j).mul( C(i,i),C(k,j) );
	tz.mul( C(k,i), C(i,j) );
	C(k,j).sub(C(k,j),tz);

	mpz_divexact(C(k,j).GetData(),C(k,j).GetData(),delta[i].GetData());
	  
      }

      C(k,i)=0;
      
    }

    delta[i+1]=C(i,i);
    
  }

   
   
  return 0;

}

int outmat(mpq_t p00,mpq_t p01,mpq_t p10,mpq_t p11) {
  
  char str[12000];
  cout << "Matrix([[ ";
  mpq_get_str (str, 10, p00);
  cout << str;
  cout << " , ";
  mpq_get_str (str, 10, p01);
  cout << str;
  cout << " ],[";
  mpq_get_str (str, 10, p10);
  cout << str;
  cout << " , ";
  mpq_get_str (str, 10, p11);
  cout << str;
  cout << " ]]);" << endl;

  return 0;
}

// **************************************************************
//
//  C, Gram, is square 
//
// updates the fraction free orthogonalization using a 2x2 matrix
// corresponding to c++ indices kappa and kappa+1 
//
// above diagonal 2 columns, diagonal block, left 2 rows 
//
// *************************************************************

int updateff(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t>& U, vector< Z_NR<mpz_t> >& delta, int kappa, ZZ_mat<mpz_t> V) {

 
  // Special case mu = 0 and V00 = 0 (sinon bug)
  //  simply a column permutation 
  // -------------------------------------------
  
  if ((V(0,0)==0) && (C(kappa,kappa+1) ==0)) {

    int i,j;

    int d;
    d=C.getCols();
    
    Z_NR<mpz_t> newd,tz;

    newd.mul(C(kappa+1,kappa+1),delta[kappa]);
    mpz_divexact(newd.GetData(),newd.GetData(),delta[kappa+1].GetData());

    for (i=0; i<d; i++) {
      tz=U(i,kappa+1);
      U(i,kappa+1)=U(i,kappa);
      U(i,kappa)=tz;
    }

    for (i=0; i<kappa; i++) {
      tz=C(i,kappa+1);
      C(i,kappa+1)=C(i,kappa);
      C(i,kappa)=tz;
    }

    for (j=kappa+2;j<d; j++) {

      tz=C(kappa,j);

      C(kappa,j).mul(C(kappa+1,j),delta[kappa]);
      mpz_divexact( C(kappa,j).GetData(), C(kappa,j).GetData(),delta[kappa+1].GetData());

      C(kappa+1,j).mul(tz,newd);
      mpz_divexact( C(kappa+1,j).GetData(), C(kappa+1,j).GetData(),delta[kappa].GetData());

    }

    C(kappa,kappa)=newd;
    delta[kappa+1]=newd;

    
  }
  else { 
   
    // Computation of the left fraction free 2x2 transformation
    // --------------------------------------------------------

    mpq_t t1,t2,t3,t4;
    mpq_init(t1);
    mpq_init(t2);
    mpq_init(t3);
    mpq_init(t4);
  
    mpq_t p00,p01,p10,p11;
    mpq_init(p00);
    mpq_init(p01);
    mpq_init(p10);
    mpq_init(p11);

 
  
    mpq_set_ui(p10,0,1);
    mpq_set_z(p01,V(1,0).GetData());

    mpq_set_z(t1,C(kappa,kappa+1).GetData());
    mpq_set_z(t2,C(kappa,kappa).GetData());
    mpq_div(p00,t1,t2); // mu up to the sign
    
    mpq_set_z(t1,V(1,0).GetData());
    mpq_set_z(t2,V(0,0).GetData());

    mpq_mul(p00,p00,t1);
    mpq_add(p00,p00,t2);
    mpq_canonicalize(p00);

    mpq_neg(p11,p00);

    mpq_inv(p11,p11);
    
  
  
    Z_NR<mpz_t> z1,z2;

    mpq_t g1,g2;
    mpq_init(g1);
    mpq_init(g2);
  
    z1.mul(C(kappa,kappa),V(0,0));
    z2.mul(C(kappa,kappa+1),V(1,0));
    z1.add(z1,z2);
    mpq_set_z(t1,z1.GetData());
    mpq_set_z(t2,delta[kappa].GetData());
    
    mpq_div(g1,t1,t2);
   
    mpq_canonicalize(g1);

    z1.mul(C(kappa+1,kappa),V(0,0));
    z2.mul(C(kappa+1,kappa+1),V(1,0));
    z1.add(z1,z2);
    mpq_set_z(t1,z1.GetData());
    mpq_set_z(t2,delta[kappa+1].GetData());
    mpq_div(g2,t1,t2);
   
    mpq_canonicalize(g2);

  
   
    mpq_mul(t1,p00,g1);
    mpq_mul(t2,p01,g2);
    mpq_add(t3,t1,t2);
    mpq_canonicalize(t3);

    mpq_mul(t1,p10,g1);
    mpq_mul(t2,p11,g2);
    mpq_add(t4,t1,t2);
    mpq_canonicalize(t4);

    mpq_set(g1,t3);
    mpq_set(g2,t4);

    mpq_div(t1,g2,g1);
    
    mpq_mul(t2,p00,t1);
    mpq_sub(p10,p10,t2);
    mpq_mul(t2,p01,t1);
    mpq_sub(p11,p11,t2);

  
  
    mpq_set_z(t1,delta[kappa].GetData());
    mpq_set_z(t2,delta[kappa+1].GetData());
    mpq_div(t1,t1,t2);
   
    mpq_mul(p01,p01,t1);
    mpq_canonicalize(p01);

    // New det
    mpq_t ndet;
    mpq_init(ndet);

    z1.mul(C(kappa,kappa),V(0,0));
    z2.mul(C(kappa,kappa+1),V(1,0));
    z1.add(z1,z2);
    mpq_set_z(g1,z1.GetData());
    mpq_canonicalize(g1);

    z1.mul(C(kappa+1,kappa),V(0,0));
    z2.mul(C(kappa+1,kappa+1),V(1,0));
    z1.add(z1,z2);
    mpq_set_z(g2,z1.GetData());
    mpq_canonicalize(g2);

  
  
    mpq_mul(t1,g1,p00);
  
    mpq_mul(t2,g2,p01);
    mpq_add(ndet,t1,t2);
    mpq_canonicalize(ndet);


    mpq_set_z(t1,delta[kappa].GetData());
    mpq_div(p10,p10,t1);
    mpq_mul(p10,p10,ndet);
    mpq_canonicalize(p10);

    mpq_set_z(t1,delta[kappa+1].GetData());
    mpq_div(p11,p11,t1);
    mpq_mul(p11,p11,ndet);
    mpq_canonicalize(p11);
  

  
  
    // The left 2x2 transformation numerators, denominator delta[kappa+1]
    // ------------------------------------------------------------------
  
    ZZ_mat<mpz_t> PN;
    PN.resize(2,2);

    mpq_set_z(t2,delta[kappa+1].GetData());
  
    mpq_set(t1,p00);
    mpq_mul(p00,p00,t2);
    mpq_canonicalize(p00);
    mpq_get_num(PN(0,0).GetData(),p00);
  
    mpq_set(t1,p10);
    mpq_mul(p10,p10,t2);
    mpq_canonicalize(p10);
    mpq_get_num(PN(1,0).GetData(),p10);

    mpq_set(t1,p01);
    mpq_mul(p01,p01,t2);
    mpq_canonicalize(p01);
    mpq_get_num(PN(0,1).GetData(),p01);

    mpq_set(t1,p11);
    mpq_mul(p11,p11,t2);
    mpq_canonicalize(p11);
    mpq_get_num(PN(1,1).GetData(),p11);


    // Update of the modified determinant
    //   and we apply the transformations on integers than exact division by den
    // -------------------------------------------------------------------
  
    Z_NR<mpz_t> den;
    den=delta[kappa+1];
  
    mpq_get_num((delta[kappa+1]).GetData(),ndet);


    int i,j;
  
    // V on the columns of C 
    // ---------------------

    ZZ_mat<mpz_t> tM;
    tM.resize(kappa+2,2);

    for (i=0; i<kappa+2; i++)
      for (j=0; j<2; j++) 
	tM(i,j)=C(i,j+kappa);
  
    matprod_in(tM,V);
  
    for (i=0; i<kappa+2; i++)
      for (j=0; j<2; j++) 
	C(i,j+kappa)=tM(i,j);

  
    
    // PN on the rows
    // --------------

    int d;
    d=C.getCols();
  
    tM.resize(2,d-kappa);

    for (i=0; i<2; i++)
      for (j=0; j<d-kappa; j++) 
	tM(i,j)=C(i+kappa,j+kappa);
  
    matprod(tM,PN,tM);

  
  
    // And exact division
    // ------------------

    for (i=0; i<2; i++)
      for (j=0; j<d-kappa; j++) 
	mpz_divexact(tM(i,j).GetData(),tM(i,j).GetData(),den.GetData());

    for (i=0; i<2; i++)
      for (j=0; j<d-kappa; j++) 
	C(i+kappa,j+kappa)=tM(i,j);
  

 
    // V on the columns of U and positivity of the diagonal of C 
    // ---------------------------------------------------------

  
    tM.resize(d,2);

    for (i=0; i<d; i++)
      for (j=0; j<2; j++) 
	tM(i,j)=U(i,j+kappa);
  
    matprod_in(tM,V);
  
    for (i=0; i<d; i++)
      for (j=0; j<2; j++) 
	U(i,j+kappa)=tM(i,j);

    if (C(kappa,kappa) < 0)
      for (j=kappa; j<d; j++)
	C(kappa,j).neg(C(kappa,j));

    if (C(kappa+1,kappa+1) < 0)
      for (j=kappa+1; j<d; j++)
	C(kappa+1,j).neg(C(kappa+1,j));

    delta[kappa+1].abs(delta[kappa+1]);

    
    // mpq clear
    // ---------
    mpq_clear(t1);
    mpq_clear(t2);
    mpq_clear(t3);
    mpq_clear(t4);
  
    mpq_clear(p00);
    mpq_clear(p01);
    mpq_clear(p10);
    mpq_clear(p11);

    mpq_clear(g1);
    mpq_clear(g2);

    mpq_clear(ndet);

  } // else of special case
  
  return 0;
   
} 


// *************************************************************
//
//  Schonhage 2x2 on the triangular gram matrix
//     diagonal is positive 
//
//    last two parameters : for stats only
//
// *************************************************************

int schon(ZZ_mat<mpz_t>& V, int& swapsloc, Z_NR<mpz_t> d1, Z_NR<mpz_t> d2,
	  Z_NR<mpz_t> B1, Z_NR<mpz_t> B2, Z_NR<mpz_t> mu, int kappa, double delta_lll,
	  int nech, vector<int>& ech) {

  setId(V);

  swapsloc=0;
  
  mpz_t mpq;
  mpz_init(mpq);

  mpz_t mpr;
  mpz_init(mpr);
  
  Z_NR<mpz_t> q,B,tz,tz1,tz2,tz3;

  Z_NR<mpz_t> one, two;
  one = 1;
  two = 2;
     
  bool stop = false;

  int nd,dd;

  nd = ((int) (delta_lll*1000000));
  dd = 1000000;

  // Main loop
  // ---------

  while (stop == false)
    {
    
    // size reduction
    // --------------
    
    if ( (mpz_sgn(mu.GetData()) * mpz_sgn(B1.GetData()) ) == 1) {
      // all positive in entry by assmumption 
      
      mpz_cdiv_qr (mpq, mpr, mu.GetData(), B1.GetData());

      mpz_set(mu.GetData(),mpr);
      mpz_set(q.GetData(),mpq);

      tz.abs(mu);
      tz.mul(tz,two);
      
      if (tz > B1) {
	mu.add(mu,B1);
	q.sub(q,one);
      }

      V(0,1).submul(V(0,0),q);
      V(1,1).submul(V(1,0),q);
          
    }
    else if ( (mpz_sgn(mu.GetData()) * mpz_sgn(B1.GetData()) ) == (-1) ) {
      // mu negative in entry by assumption 
      
      mpz_fdiv_qr (mpq, mpr, mu.GetData(), B1.GetData());

      mpz_set(mu.GetData(),mpr);
      mpz_set(q.GetData(),mpq);
      
      tz.abs(mu);
      tz.mul(tz,two);
       
      if (tz > B1) {
	mu.sub(mu,B1);
	q.add(q,one);
      }

      V(0,1).submul(V(0,0),q);
      V(1,1).submul(V(1,0),q);
    }

    
    // Lovasz test
    // -----------

    // Swap

    tz1.mul(B1,B1);
    
    tz2.mul(B2,d1);
    tz3.mul(mu,mu);
    tz2.add(tz2,tz3);

    tz1.mul_ui(tz1,nd);
    tz2.mul_ui(tz2,dd);

   
    if (tz1 > tz2) {

      swapsloc+=1;
      
      B.mul(mu,mu);
      tz.mul(B2,d1);
      B.add(tz,B);
      mpz_divexact(B.GetData(),B.GetData(),B1.GetData());


      // STATS
      // -----

      {
	mpq_t qq,rho;
	mpq_init(qq);
	mpq_init(rho);

	mpq_set_z(rho,B.GetData());
	mpq_set_z(qq,B1.GetData());

	mpq_div(rho,rho,qq); // mu up to the sign

	double tf;
	tf=mpq_get_d(rho);

       	
	ech[((int) (tf*nech))]+=1; 
	  
	  //cout << "******* Ratio dyn at " << kappa << " : " << mpq_get_d(rho) << endl; 

	mpq_clear(qq);
	mpq_clear(rho);
	
      } // End stats 
      
      B1=B;
      d2=B;

      tz=V(0,1);
      V(0,1)=V(0,0);
      V(0,0)=tz;

      tz=V(1,1);
      V(1,1)=V(1,0);
      V(1,0)=tz;

     
      
    }
    else {
      stop = true; 
    }
    
  } // 2x2 LLL end loop  

 
  
  mpz_clear(mpq);
  mpz_clear(mpr);
  
  return(0);
  
}


// *************************************************************
//
//  Parameters are allocated outside
//
//  Size reduction for the two modified columns kappa and kappa+1
//
// *************************************************************



int ffreduce(ZZ_mat<mpz_t>& C,  ZZ_mat<mpz_t>& U, int kappa) {

    
  int d;
  d=C.getCols();

  int ck,i,k;

  mpz_t mpq;
  mpz_init(mpq);

  mpz_t mpr;
  mpz_init(mpr);
  
  Z_NR<mpz_t> q,r; 

  Z_NR<mpz_t> one, two;
  one = 1;
  two = 2;

  Z_NR<mpz_t> mu, diag,tz;
  
  // Loop on the two columns to reduce
  // ---------------------------------
  
  for (ck=kappa; ck<kappa+2; ck++) {

    // loop within the column 
    for (i=ck-1; i>=0; i--) {

      // compare C(i,ck) with C(i,i)

      mu = C(i,ck);
      diag =  C(i,i);
      
      // quotient 
      // --------
    
      if ( (mpz_sgn(mu.GetData()) * mpz_sgn(diag.GetData()) ) == 1) {
	// all positive in entry by assmumption 
	
	mpz_cdiv_qr (mpq, mpr, mu.GetData(), diag.GetData());

	mpz_set(mu.GetData(),mpr);
	mpz_set(q.GetData(),mpq);

	tz.abs(mu);
	tz.mul(tz,two);
      
	if (tz > diag) {
	  mu.add(mu,diag);
	  q.sub(q,one);
	}


	// Not use variable index i here !!
	// Application of the transformation on C and U 
	
	for (k=i; k>=0; k--)
	  C(k,ck).submul(C(k,i),q);
	
	for (k=0; k<d; k++)
	  U(k,ck).submul(U(k,i),q);

      }
      else if ( (mpz_sgn(mu.GetData()) * mpz_sgn(diag.GetData()) ) == (-1) ) {
	// mu negative in entry by assumption 
      
	mpz_fdiv_qr (mpq, mpr, mu.GetData(), diag.GetData());
	
	mpz_set(mu.GetData(),mpr);
	mpz_set(q.GetData(),mpq);
	
	tz.abs(mu);
	tz.mul(tz,two);
	
	if (tz > diag) {
	  mu.sub(mu,diag);
	  q.add(q,one);
	}

	// Not use variable index i here !!
	// Application of the transformation on C and U 

	for (k=i; k>=0; k--)
	  C(k,ck).submul(C(k,i),q);
	
	for (k=0; k<d; k++)
	  U(k,ck).submul(U(k,i),q);
	 
      }

    } // end within the column 
    

  } // end loop 2 cols 

  
 
  return(0); 
}




// *************************************************************
//
// LLL reduction : computation of the transformation matrix U
//
// *************************************************************

int fflll(ZZ_mat<mpz_t>& U,  ZZ_mat<mpz_t> A, double delta_lll) {

  // for STATS
  int nech=100;
  vector< int > ech;
  ech.resize(nech);

  for (int i=0; i<nech; i++)
    ech[i]=0;
  
  int d;

  setId(U);
  
  d=A.getCols();
    
  ZZ_mat<mpz_t> C;
  C.resize(d,d);

  vector< Z_NR<mpz_t> > delta;
  delta.resize(d);

  ff(C,delta,A);

  // For Heckler & Thiele (log_2, c=2)
  // ---------------------------------
  
  vector<FP_NR<mpfr_t> > v;
  v.resize(d);

  vector<FP_NR<mpfr_t> > cv; // Constant part of v_i
  cv.resize(d);

  Z_NR<mpz_t> tz;
  FP_NR<mpfr_t> tf,tf2;
  
  for (int i=1; i<=d; i++) {
   
    tz=-i*(d-i);
    tf.set_z(tz);
    tf.mul_2si(tf,-1);

    cv[i-1]=tf;
    
    tz=C(d-1,d-1);
    tf.set_z(tz);
    mpfr_log2(tf.get_data(),tf.get_data(),GMP_RNDN); 
 
    tz=-i;
    tf2.set_z(tz);
    tf.mul(tf,tf2);
    
    tz=d;
    tf2.set_z(tz);

    tf.div(tf,tf2);

    cv[i-1].add(cv[i-1],tf);

  }
  // --------- end HT

  
  
  // Unimodular 2x2 transformation
  ZZ_mat<mpz_t> V;
  V.resize(2,2);

  // ------------------------------------
  // Main dynamical LLL reduction loop
  // ------------------------------------

  int kappa;

  int K;

  int nbswaps,swapsiter,swapsloc;

  nbswaps=0;
  swapsiter=0;

  int bounditer;

  bounditer = d*d;
  
  // Main LLL loop
  // -------------

  for (K=0; K<bounditer; K++) {

    // Loop on the 2x2 Gauss application
    // ---------------------------------
    
    for (kappa=d-2; kappa >=0; kappa--) {

      schon(V, swapsloc, delta[kappa], delta[kappa+1], C(kappa,kappa),C(kappa+1,kappa+1),C(kappa,kappa+1),
	    kappa,delta_lll,nech,ech);

          
      nbswaps += swapsloc;
      swapsiter += swapsloc;

     
      
      updateff(C,U,delta,kappa,V);

      ffreduce(C,U,kappa);

    }

    cout << "Iteration " << K << endl;
    cout << "       swaps : " << swapsiter << endl;
    cout << "       total : " << nbswaps << endl;

    
    
    // Heckler & Thiele
    // ----------------
    for (int i=0; i<d; i++) {

      tz=C(i,i);
      tz.abs(tz);
      tf.set_z(tz);
      mpfr_log2(tf.get_data(),tf.get_data(),GMP_RNDN); 

      v[i]=tf;
      v[i].add(cv[i],v[i]);

    
    }

    tf=v[0];
    
    for (int i=1; i<d; i++)
      if (tf < v[i]) tf=v[i];
    mpfr_exp2(tf.get_data(),tf.get_data(),GMP_RNDN); 

    cout  << "       HT " << tf << endl << endl;
    
    // ------- end HT

    // Early termination
    // -----------------
    if (swapsiter==0) {
      cout << endl; 
      cout << "Early termination after " << K+1 << " iterations" << endl;
      K=bounditer; 
    }

    swapsiter=0;


  }  // end main LLL loop

  // STATS
  cout << endl << endl << "Dyn repartition " << endl << endl;

  cout << "[";
  for (int i=0; i<nech-2; i++)
    cout << "[ " << i << "," << ech[i] << " ],";
  cout << "[ " << nech-2 << "," << ech[nech-2] << " ]";
  cout << "]" << endl; 
  
  // Final size reduction ??
  // --------------------

  //print2maple(U,d,d);
  

  return nbswaps;
  
} 

// ***********************************************************
       
int main(int argc, char *argv[]) {
  
  ZZ_mat<mpz_t> A;

  int n,d;
  double delta_lll;

  command_line_basis(A, n, d, delta_lll, argc, argv); 

  
  ZZ_mat<mpz_t> U;
  U.resize(d,d);

  int dynswaps;
  
  dynswaps=fflll(U,A,delta_lll);

  
  Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> > B(A,NO_TRANSFORM,DEF_REDUCTION);

  B.hlll(delta_lll);

  // Stats
  // -----
  
  double ta,ua,va,wa;
  
  matprod_in(A,U);   // *** erases A
  hplll::ratio<mpz_t>(A,ta,ua,va,wa);
  
  double tb,ub,vb,wb;
  
  hplll::ratio<mpz_t>(B.getbase(),tb,ub,vb,wb);

  cout << endl << ".. log 2 Frobenius norm cond: " << ta << endl;
  cout         << "                              " << tb << endl;
  
  cout << ".. Average diagonal ratio: " << ua << endl;
  cout << "                           " << ub << endl;
  
  cout << ".. Max diagonal ratio: " << va << endl;
  cout << "                       " << vb << endl;
  
  cout << ".. First vector quality: " << wa << endl;
  cout << "                         " << wb << endl;

  
  cout << endl << "Dyn swaps: " << dynswaps << endl;
  cout << endl << "LLL swaps: " << B.nbswaps << endl << endl;

  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A,NO_TRANSFORM,DEF_REDUCTION);
  verboseDepth=0;
  T.isreduced(delta_lll-0.1);
  
}
