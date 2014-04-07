// Copyright (C) IBM, All Rights Reserved

/*** fhe.h - header for implementation of (a variant of) Gentry's
 *           Fully-homomorphic scheme
 *
 *   - FHEkeys::keyGen():  generate public/secret keys
 *   - FHEkeys::encrypt(): encrypt a bit
 *   - FHEkeys::decrypt(): recover the plaintext bit from the ciphertext
 *   - FHEkeys::recrypt(): perform re-crypt (homomorphic decryption) to refersh a cipehrtext
 *
 ********************************************************************/
#ifndef _FHE_H_
#define _FHE_H_
#include <vector>
#include <cmath>

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZX.h>

typedef ZZ FHEctxt;

class FHEparams {
 public:
  unsigned long secprm;// the security parameter
  double mu;           // the BDD-hardness parameter

  unsigned long s;     // sparse-subset size
  unsigned long S;     // big-set size
  unsigned long p;     // the precision parameter
  unsigned long t;     // bitsize of coeffs in generating polynomial
  unsigned long logn;  // log the dimenssion
  unsigned long logR;  // log the ratio between elements in a big set
  unsigned long noise; // expected # of nonzero entries in fresh cipehrtexts

  FHEparams(unsigned long secprm=64, unsigned long t=0, unsigned long logn=0)
    {  setPrms(secprm, t, logn);  }

  long operator==(const FHEparams& other) const
    { return (secprm == other.secprm
	      && mu == other.mu
	      && s == other.s
	      && S == other.S
	      && p == other.p
	      && t == other.t
	      && logn == other.logn
	      && logR == other.logR
	      && noise == other.noise);}
  // These functions are defined in fhe-utils.cc
  void setPrms(unsigned long sec, unsigned long t=0, unsigned long logn=0);
  void print(std::ostream& st) const;
};

/* The bulk of the "squashed" public key consists of s blocks, where s
 * is the size of the small subset. Each block is logically a sequence
 * of S elements xi modulo det (where S is the "big-set size") and
 * S ciphertexts, encrypting S of the secret key bits.
 *
 * To save on space, we only keep one elements x for every block and
 * sets xi = R^i * x mod det, for the parameter R.  Also, we only keep
 * about sqrt(2S) ciphertexts, two of which are encryptions of "1" and
 * the others encryptions of "0", such that the encryption of every
 * secret-key bit is a product of two of them.
 *
 * Within each block we have exactly one xi that participates in the
 * sparse subset, and this xi is indexed by the pair of 1-ciphertexts.
 *
 * For example, if S=2^13 then we have 128 ciphertexts for the block,
 * and if (i1,i2) are the indexes of the two that encrypt a "1", then
 * the element from that block that participates in the sparse subset
 * is xi, where i is
 *    i = (i1-1)*128 -(i1 choose 2) +(i2-i1) = i1*(255-i1)/2 +i2 -128
 */
class PKblock {
 public:
  ZZ x;         // the sequence of elements is xi = x*2^i mod det
  long idx;     // index of the xi that belongs to the sparse subset

  long operator==(const PKblock& other) const
    { return (x == other.x && idx == other.idx); }
};

class FHEkeys {
 public:
  FHEparams prms;   // Parameters
  ZZ det, root, w;  // Public/secret key for underlying "somehwat" scheme
  std::vector<PKblock> pkBlocks;  // The squeash-enabling parts
  vec_ZZ ctxts;     // "compact" encryption of secret-key bits

  FHEkeys(unsigned long sec=64, unsigned long t=0, unsigned long logn=0):
    prms(sec,t,logn) {}

  // Utility methods

  const FHEparams& getPrms() const {return prms;}
  void setPrms(unsigned long sec, unsigned long t=0, unsigned long logn=0)
    { prms.setPrms(sec,t,logn); }

  const ZZ& getDet() const {return det;}
  const ZZ& getRoot() const {return root;}
  const ZZ& getSkey() const {return w;}

  long operator==(const FHEkeys& other) const
    { return (prms == other.prms
	      && det == other.det
	      && root== other.root
	      && w   == other.w
	      && pkBlocks == other.pkBlocks
	      && ctxts == other.ctxts);     }
  
  // The next few functions are defined in fhe-utils.cc
  friend std::istream& operator>>(std::istream& s, FHEkeys& keys);
  friend std::ostream& operator<<(std::ostream& s, FHEkeys& keys);
  void inputKeys(std::istream& s, int inIdxs=1, int inW=1);
  void outputKeys(std::ostream& s, int outIdxs=1, int outW=1);

  // returns (an approx of) log(det/(c*w)). This should
  // be at least p+1 if we want to be able to reCrypt c.
  long testCtxt(const FHEctxt& c) const;

  // returns in j1,j2 the i'th pair (in lexicographic order) in {m choose 2}
  static void encodeIndex(unsigned long& j1, unsigned long& j2,
			  unsigned long i, unsigned long m);

  // returns the smallest integer m s.t. m(m-1)/2 >= S
  static unsigned long mChoose2(unsigned long S) {
    return (unsigned long) ceil(2*sqrt((double)S));
  }

  // Real work

  int keyGen();   // returns # of trials until success (in fhe-keygen.cc)

  // return TRUE on success, FALSE otherwise (defined in fhe-enc.cc)
  bool encrypt(vec_ZZ& c, unsigned int m[], int num) const; // batch encryption
  bool encrypt(FHEctxt& c, unsigned int m) const {
    vec_ZZ vc(INIT_SIZE, 1);
    vc[0] = c;
    if (!encrypt(vc, &m, 1)) return false;
    c = vc[0];
    return true;
  }

  // decrypt a single bit, returns the plaintext bit or -1 on error
  // (defined in fhe-dec.cc)
  unsigned int decrypt(const FHEctxt& c) const; 

  // just for symmetry, a trivial batch decryption,
  // returns the number of bits that were decrypted successfully
  int decrypt(unsigned int m[], const vec_ZZ& c) const {
    int i;
    for (i=0; i<c.length(); i++) {
      if ((m[i]=decrypt(c[i])) ==-1) break;
    }
    return i; // returns the number of successfully decrypted bits
  }

  void recrypt(FHEctxt& c) const;
 private:
  // called from recrypt to process one public-key block
  void processBlock(vec_ZZ& vars, const FHEctxt& c, long i) const; 
};
#endif /* ifndef _FHE_H_ */
