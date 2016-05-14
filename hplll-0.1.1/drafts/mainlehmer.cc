/* Integer matrix nullspace test file  

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

  typedef __int128_t ZT;

  Z_NR<mpz_t> a,b;

  a=25443;
  a.mul(a,a);
  a.mul(a,a);
  a.mul(a,a);   // 175608908455044406357860029088740001   size: 117

  cout << "a: " << a << "   size: " << size_in_bits(a) << endl;

  b=232;
  b.mul(b,b);
  b.mul(b,b);
  b.mul(b,b);
  b.mul(b,b);
  b.neg(b);     // -70438120351099559671412028074440523776   size: 125
  
  cout << "b: " << b << "   size: " << size_in_bits(b) << endl;

  // MPZ to 128
  
  char str[sizeof(__int128_t)+4];

  mpz_get_str (str, 10, b.getData());

  cout << "str " << str << endl;

  std::string str_dec = str;
  
  //std::string str_hex = "40c3";
  //std::string str_bin = "-10010110001";
  //std::string str_auto = "0x7f";

  std::string::size_type sz;   // alias of size_t

  __int128_t c;
  
  c = std::stoi (str,&sz);


 //  int print_uint128(uint128_t n) {
//   if (n == 0)  return printf("0\n");

//   char str[40] = {0}; log10(1 << 128) + '\0'
//   char *s = str + sizeof(str) - 1; start at the end
//   while (n != 0) {
//     if (s == str) return -1; never happens

//     *--s = "0123456789"[n % 10]; save last digit
//     n /= 10;                     drop it
//   }
//   return printf("%s\n", s);
// }
  
  // 128 to mpz via  mpz_set_str ds nr-ld 
  
  return 0;
}
