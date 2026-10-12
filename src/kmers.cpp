#include <stdlib.h>
#include "dada.h"
#include "kmers_simd.h"
// #if _WIN32
// #include <intrin.h>
// #endif
// [[Rcpp::interfaces(cpp)]]

/************* KMERS *****************
 * Current Kmer implementation assumes A/C/G/T only.
 * Current Kmer implementation requires kmer size of 3-7 only.
 */

// DadaS: the implementations live in kmers_simd.h (SSE2 on x86, NEON on ARM, scalar elsewhere).
// These wrappers keep the original interface.
double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  return kmer_dist_scalar(kv1, len1, kv2, len2, k);
}

double kord_dist(uint16_t *kord1, int len1, uint16_t *kord2, int len2, int k) {
  return kord_dist_scalar(kord1, len1, kord2, len2, k);
}

#if DADAS_SIMD
double kmer_dist_SSEi(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  return kmer_dist_simd16(kv1, len1, kv2, len2, k);
}
// Uses kmers packed into 8 bit integers. Can overflow, in which case negative value is returned.
double kmer_dist_SSEi_8(uint8_t *kv1, int len1, uint8_t *kv2, int len2, int k) {
  return kmer_dist_simd8(kv1, len1, kv2, len2, k);
}
double kord_dist_SSEi(uint16_t *kord1, int len1, uint16_t *kord2, int len2, int k) {
  return kord_dist_simd(kord1, len1, kord2, len2, k);
}
#else // no SIMD: never called (Rmain.cpp sets SSE=0), defined for completeness
double kmer_dist_SSEi(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  return kmer_dist_scalar(kv1, len1, kv2, len2, k);
}
double kmer_dist_SSEi_8(uint8_t *kv1, int len1, uint8_t *kv2, int len2, int k) {
  return -1.; // forces the 16-bit version
}
double kord_dist_SSEi(uint16_t *kord1, int len1, uint16_t *kord2, int len2, int k) {
  return kord_dist_scalar(kord1, len1, kord2, len2, k);
}
#endif

void assign_kmer8(uint8_t *kvec8, const char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  size_t len = strlen(seq);
  if(len <= 0 || len > SEQLEN) { dada_stop("Unexpected sequence length."); }
  if(k >= len || k < 3 || k > 8) { dada_stop("Invalid kmer-size."); }
  size_t klen = len - k + 1; // The number of kmers in this sequence
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  uint16_t *kvec = (uint16_t *) malloc(n_kmers * sizeof(uint16_t)); //E
  if (kvec == NULL)  dada_stop("Memory allocation failed.");
  for(kmer=0;kmer<n_kmers;kmer++) { kvec[kmer] = 0; }
  
  if(len <=0 || len > SEQLEN) {
    dada_stop("Unexpected sequence length.");
  }
  
  for(i=0; i<klen; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        dada_stop("Unexpected nucleotide.");
        kmer = 999999;
        break;
      }
      kmer = 4*kmer + nti;
    }
    
    // Make sure kmer index is valid. This doesn't solve the N's/-'s
    // issue though, as the "length" of the string (# of kmers) needs
    // to also reflect the reduction from the N's/-'s
    if(kmer == 999999) { ; } 
    else if(kmer >= n_kmers) {
      dada_stop("Kmer index out of range.");
    } else { // Valid kmer
      kvec[kmer]++;
    }
  }
  // Move to destination uint8_t vector w/ overflow checks
  for(kmer=0;kmer<n_kmers;kmer++) { 
    if(kvec[kmer]<255) {
      kvec8[kmer] = (uint8_t)kvec[kmer];
    } else {
      kvec8[kmer] = 255;
    }
  }
  free(kvec);
}

void assign_kmer(uint16_t *kvec, const char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  size_t len = strlen(seq);
  if(len <= 0 || len > SEQLEN) { dada_stop("Unexpected sequence length."); }
  if(k >= len || k < 3 || k > 8) { dada_stop("Invalid kmer-size."); }
  size_t klen = len - k + 1; // The number of kmers in this sequence
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  for(kmer=0;kmer<n_kmers;kmer++) { kvec[kmer] = 0; }
  
  if(len <=0 || len > SEQLEN) {
    dada_stop("Unexpected sequence length.");
  }
  
  for(i=0; i<klen; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        dada_stop("Unexpected nucleotide.");
        kmer = 999999;
        break;
      }
      kmer = 4*kmer + nti;
    }
    
    // Make sure kmer index is valid. This doesn't solve the N's/-'s
    // issue though, as the "length" of the string (# of kmers) needs
    // to also reflect the reduction from the N's/-'s
    if(kmer == 999999) { ; } 
    else if(kmer >= n_kmers) {
      dada_stop("Kmer index out of range.");
    } else { // Valid kmer
      kvec[kmer]++;
    }
  }
}

// Returns the vector of the kmers in order
void assign_kmer_order(uint16_t *kord, char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  size_t len = strlen(seq);
  if(len <=0 || len > SEQLEN) { dada_stop("Unexpected sequence length."); }
  if(k >= len || k < 1 || k > 8) { dada_stop("Invalid kmer-size."); }
  size_t klen = len - k + 1; // The number of kmers in this sequence
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  if (kord == NULL)  dada_stop("Memory allocation failed.");
  for(i=0;i<klen;i++) { kord[i] = 0; }
  
  for(i=0; i<klen; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        dada_stop("Unexpected nucleotide.");
        kmer = 999999;
        break;
      }
      kmer = 4*kmer + nti;
    }
    
    // Make sure kmer index is valid. This doesn't solve the N's/-'s
    // issue though, as the "length" of the string (# of kmers) needs
    // to also reflect the reduction from the N's/-'s
    if(kmer == 999999) { ; } 
    else if(kmer >= n_kmers) { // maybe check if above uin16_t max? Defined in header <stdint.h>, UINT16_MAX
      dada_stop("Kmer index out of range.");
    } else { // Valid kmer
      kord[i] = kmer;
    }
  }
}

