
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
  
  return 0;
   
} 


// *************************************************************
//
//  Schonhage 2x2 on the triangular gram matrix
//     diagonal is positive 
//
// *************************************************************

int schon(ZZ_mat<mpz_t>& V, Z_NR<mpz_t> d1, Z_NR<mpz_t> d2, Z_NR<mpz_t> B1, Z_NR<mpz_t> B2, Z_NR<mpz_t> mu, double delta_lll) {

  setId(V);
  
  mpz_t mpq;
  mpz_init(mpq);

  mpz_t mpr;
  mpz_init(mpr);
  
  Z_NR<mpz_t> q,B,tz,tz1,tz2,tz3;

  Z_NR<mpz_t> one, two;
  one = 1;
  two = 2;
     
  bool stop = false;
  
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

    }

    
    V(0,1).submul(V(0,0),q);
    V(1,1).submul(V(1,0),q);

    

    // Lovasz test
    // -----------

    // Swap

    tz1.mul(B1,B1);
    
    tz2.mul(B2,d1);
    tz3.mul(mu,mu);
    tz2.add(tz2,tz3);
   
    
    if (tz1 > tz2) {

      B.mul(mu,mu);
      tz.mul(B2,d1);
      B.add(tz,B);
      mpz_divexact(B.GetData(),B.GetData(),B1.GetData());

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



int ffreduce(ZZ_mat<mpz_t>& U,  ZZ_mat<mpz_t> C, int kappa) {

  int d;
  d=C.getCols();

  int ck,i;

  mpz_t mpq;
  mpz_init(mpq);

  mpz_t mpr;
  mpz_init(mpr);
  
  Z_NR<mpz_t> q,r; 

  Z_NR<mpz_t> one, two;
  one = 1;
  two = 2;

  Z_NR<mpz_t> mu, diag;
  
  // Loop on the two columns to reduce
  // ---------------------------------
  
  for (ck=kappa; ck<kappa+2; col++) {

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
	
      }






      
      

    } // end within the column 
    

  } // end loop 2 cols 
  
    <// ET update de la transformation 
}




// *************************************************************
//
// LLL reduction : computation of the transformation matrix U
//
// *************************************************************

int fflll(ZZ_mat<mpz_t>& U,  ZZ_mat<mpz_t> A, double delta_lll) {

  int d;

  setId(U);
  
  d=A.getCols();
    
  ZZ_mat<mpz_t> C;
  C.resize(d,d);

  vector< Z_NR<mpz_t> > delta;
  delta.resize(d);

  ff(C,delta,A);
  
  // Unimodular 2x2 transformation
  ZZ_mat<mpz_t> V;
  V.resize(2,2);

  // ------------------------------------
  // Main dynamical LLL reduction loop
  // ------------------------------------

  int kappa;

  int K;
  
  // Main LLL loop
  // -------------

  for (K=0; K<d*d; K++) {

    // Loop on the 2x2 Gauss application
    // ---------------------------------
    
    for (kappa=d-2; kappa >=0; kappa--) {
      
      schon(V, delta[kappa], delta[kappa+1], C(kappa,kappa),C(kappa+1,kappa+1),C(kappa,kappa+1),delta_lll);

      updateff(C,U,delta,kappa,V);

      print2maple(U,d,d);
      
      print2maple(C,d,d);

      // Size reduction columns kappa and kappa+1
      //  and corresponding update of U 
      // Update U, par V ** modulo ? **

    }

  }
  // end main LLL loop
  // *****************
  

  return 0;
  
} 

// ***********************************************************
       
int main(int argc, char *argv[]) {
  
  ZZ_mat<mpz_t> A;

  int n,d;
  double delta_lll;

  command_line_basis(A, n, d, delta_lll, argc, argv); 

  
  ZZ_mat<mpz_t> U;
  U.resize(d,d);

  fflll(U,A,delta_lll);
  
 
  
}
