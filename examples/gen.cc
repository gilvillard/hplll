
#include <hplll.h>

int main(int argc, char *argv[]) {

  ZZ_mat<mpz_t> A;

  int n,d;
  double delta;

  command_line_basis(A, n, d, delta, argc, argv); 

  cout << A << endl;

  
  }
