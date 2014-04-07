// Copyright (C) IBM, All Rights Reserved

#include <cstring>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <NTL/tools.h>
#include <NTL/ZZ.h>

NTL_CLIENT
#include "fhe.h"

/* Modes of call:
 *   ops = 0: only read key from input (used to test I/O)
 *         1: only generate key
 *         2: generate key and perform re-crypt test
 *
 *   mode = 0: write (read) entire key
 *          1: write (read) w, but not the indexes
 *          2: write (read) only public key, without w or indexes
 ********************************************************************/
int  main(int argc, char *argv[])
{
  cerr << "  Usage: "<<argv[0]<< " ops iomode [bitsize dimension [prg-seed]]\n\n";
  cerr << "    ops = 0: only read key from input (used to test I/O)\n";
  cerr << "          1: only generate key\n";
  cerr << "          2: generate key and perform re-crypt test\n";
  cerr << " iomode = 0: write (read) entire key\n";
  cerr << "          1: write (read) w, but not the indexes\n";
  cerr << "          2: write (read) only public key, without w or indexes\n";

  if (argc==2 && strcmp("--help",argv[1])==0) {
    exit(0);
  }

  long logDim = 11;
  long nBits = 384;
  int ops=2, mode=0;

  if (argc>1) ops =atoi(argv[1]);
  if (argc>2) mode=atoi(argv[2]);
  if (argc>3 && atol(argv[3])>1) nBits=atol(argv[3]);
  if (argc>4 && atol(argv[4])>1) logDim=NextPowerOfTwo(atol(argv[4]));

  // If user specified PRG-seed then use it, else use current time
  ZZ prgSeed;
  if (argc>5) conv(prgSeed,argv[5]);
  else conv(prgSeed,time (NULL));
  cerr << endl << "PRG initialized with " << prgSeed << endl;
  SetSeed(prgSeed);        // initilize the RNG

  FHEkeys keys(72,nBits,logDim);  // Note: we use 72-bit security level
  cerr << "------------------------------------\n";
  if (ops > 0) { // generate the key yourself
    keys.getPrms().print(cerr);
    keys.keyGen();
  }
  else {         // read key from cin
    if (mode == 0) cin >> keys;                 // read entire key
    else if (mode == 1) keys.inputKeys(cin, 0); // read w, not indexes
    else keys.inputKeys(cin, 0, 0);             // read only public key
    keys.getPrms().print(cerr);
  }
  cerr << "------------------------------------\n";

  // write indexes to cerr (this is the solution to the challenge)
  cerr << "Challenge solution: indexes are:\n";
  if (ops > 0 || mode < 2) for (int i=0; i<keys.pkBlocks.size(); i++) { 
      cerr << " " << setw(4) << keys.pkBlocks[i].idx;
    }
  cerr << endl;

  ZZ::HexOutput = 0;  // NTL modification: I/O in hexadecimal
  if (mode == 0) cout << keys;               // write entire key
  else {
    if (mode == 1) keys.outputKeys(cout, 0); // write w, not indexes
    else keys.outputKeys(cout, 0, 0);        // write only public key
  }
  cerr << "------------------------------------\n";

  if (ops<2) exit(0); // should we proceed to re-crypt test?

  // Generate a fresh encryption of one
  FHEctxt c;
  keys.encrypt(c,1);
  cerr << "c = Enc(1)\n";

  // Square the ciphertext, then apply recrypt and check the result
  // is still decrypted as 1
  for (int i=0; i<8; i++) {
    SqrMod(c,c,keys.getDet());
    if (keys.testCtxt(c) < keys.getPrms().p +1) // sanity check
      cerr << "c * w larger than det / 2^{p+1}!!\n";
    keys.recrypt(c);
    cerr << "c = reCrypt(c^2) decrypted as " << keys.decrypt(c) <<endl;
  }
}
