
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
// updates the fraction free orthogonalization using a 2x2 matrix
// corresponding to c++ indices kappa and kappa+1 
//
// above diagonal 2 columns, diagonal block, left 2 rows 
//
// *************************************************************

int updateff(ZZ_mat<mpz_t>& C, vector< Z_NR<mpz_t> >& delta, int kappa, ZZ_mat<mpz_t> U) {


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
  mpq_set_z(p01,U(1,0).GetData());

  mpq_set_z(t1,C(kappa,kappa+1).GetData());
  mpq_set_z(t2,C(kappa,kappa).GetData());
  mpq_div(p00,t1,t2); // mu up to the sign

  mpq_set_z(t1,U(1,0).GetData());
  mpq_set_z(t2,U(0,0).GetData());

  mpq_mul(p00,p00,t1);
  mpq_add(p00,p00,t2);
  mpq_canonicalize(p00);

  mpq_neg(p11,p00);
  mpq_inv(p11,p11);

  Z_NR<mpz_t> z1,z2;

  mpq_t g1,g2;
  mpq_init(g1);
  mpq_init(g2);
  
  z1.mul(C(kappa,kappa),U(0,0));
  z2.mul(C(kappa,kappa+1),U(1,0));
  z1.add(z1,z2);
  mpq_set_z(t1,z1.GetData());
  mpq_set_z(t2,delta[kappa].GetData());
  mpq_div(g1,t1,t2);
  mpq_canonicalize(g1);

  z1.mul(C(kappa+1,kappa),U(0,0));
  z2.mul(C(kappa+1,kappa+1),U(1,0));
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

  z1.mul(C(kappa,kappa),U(0,0));
  z2.mul(C(kappa,kappa+1),U(1,0));
  z1.add(z1,z2);
  mpq_set_z(g1,z1.GetData());
  mpq_canonicalize(g1);

  z1.mul(C(kappa+1,kappa),U(0,0));
  z2.mul(C(kappa+1,kappa+1),U(1,0));
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
  
  outmat(p00,p01,p10,p11);
	  
  // cas special indice 0 ?
  // faire colonnes avant la diag

  // Mettre à jour delta 


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

// ***********************************************************
       
int main(int argc, char *argv[]) {
  
  ZZ_mat<mpz_t> A;

  int n,d;
  double delta_lll;

  command_line_basis(A, n, d, delta_lll, argc, argv); 

  cout << A << endl << endl;

  ZZ_mat<mpz_t> C;
  C.resize(d,d);

  vector< Z_NR<mpz_t> > delta;
  delta.resize(d);

  ff(C,delta,A);

  // Unimodular 2x2 transformation
  ZZ_mat<mpz_t> U;
  U.resize(2,2);
  U(0,0)=3; U(0,1)=1; U(1,0)=-32; U(1,1)=-11;

  
  
  updateff(C,delta,4,U);
  
 
  
}
