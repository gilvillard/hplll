/* Sam 14 mai 2016 17:12:58 CEST

With 64 bits long int 
Conversions between mpz_t and __int128_t

Created Dim  7 avr 2013 16:54:03 CEST
Copyright (C) 2013      Gilles Villard 

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


#include "hplll.h"


using namespace hplll; 


/* ***********************************************

          MAIN   

   ********************************************** */


int main(int argc, char *argv[])  {

 
  int status =0;

  int k,i;

  for (i=1; i<1000; i++) {

    // Positive 
    for (k=1; k<=127; k++) {

      Z_NR<mpz_t> a,b,one;
      one = 1;

      Z_NR<__int128_t> r,oone;
      oone = 1;

      a.randb(k);

      mpz_set_128int(r,a); 

      r.sub(r,oone);

      mpz_get_128int(b,r);

      a.sub(a,one);

      status |= a.cmp(b); 

    }
  
    // Negative  
    for (k=1; k<=127; k++) {

      Z_NR<mpz_t> a,b,one;
      one = 1;

      Z_NR<__int128_t> r,oone;
      oone = 1;

      a.randb(k);

      a.neg(a);

      mpz_set_128int(r,a); 

      r.add(r,oone);

      mpz_get_128int(b,r);

      a.add(a,one);

      status |= a.cmp(b); 

   
    }
  } // i: nb of loop tests 

  return status;
}
