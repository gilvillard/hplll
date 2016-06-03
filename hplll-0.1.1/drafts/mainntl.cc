/* NTL tests   

Created Jeu  2 jui 2016 16:28:42 CEST  
Copyright (C) 2016      Gilles Villard 

This file is part of the hplll Library 

The hplll Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The hplll Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPO
SE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the hplll Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */


#include <sstream>
#include <iostream>

#include "hplll.h"

#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/LLL.h>

//#include "ideal.h"


/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

using namespace NTL; 


int main(int argc, char *argv[])  {
  
  Mat<ZZ> B(INIT_SIZE, 2,2) ;
  int d,dbits; 

  ZZ_mat<mpz_t> AT;  // fpLLL  

  filebuf fb;
  fb.open ("ntl.txt",ios::in);
  iostream os(&fb);
  os >>  B ;
  fb.close();

  
  cout << B << endl;
  
  
 
  return 0;
}
