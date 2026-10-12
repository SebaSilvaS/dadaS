// DadaS: distancias de k-mers con SIMD para x86 (SSE2) y ARM (NEON), más la versión escalar.
//
// Este archivo no depende de R ni de Rcpp, para poder compilarlo también en el banco de pruebas
// (workflow/bench/). Lo incluye kmers.cpp.
//
// En DADA2 las versiones vectorizadas solo existían para x86_64; en ARM (Apple Silicon, Graviton,
// servidores ARM) se usaba la versión escalar. Aquí las versiones NEON replican exactamente la
// aritmética de las SSE2 (mismos carriles, mismas sumas saturadas o con desborde, misma detección
// de desborde), así que dan resultados idénticos bit a bit.
//
// DADAS_SIMD es 1 si hay una implementación vectorizada para esta arquitectura.
#ifndef DADAS_KMERS_SIMD_H
#define DADAS_KMERS_SIMD_H

#include <stdint.h>
#include <stddef.h>

#if defined(__x86_64__) || defined(__x86_64)
#include <emmintrin.h>
#define DADAS_SIMD 1
#elif defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#define DADAS_SIMD 1
#else
#define DADAS_SIMD 0
#endif

/* ---------- Escalar (referencia, cualquier arquitectura) ---------- */

static inline double kmer_dist_scalar(const uint16_t *kv1, int len1, const uint16_t *kv2, int len2, int k) {
  int n_kmer = 1 << (2*k); // 4^k kmers
  uint16_t dotsum = 0;
  for(int i=0; i<n_kmer; i++) {
    dotsum += (kv1[i] < kv2[i] ? kv1[i] : kv2[i]);
  }
  double dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}

static inline double kord_dist_scalar(const uint16_t *kord1, int len1, const uint16_t *kord2, int len2, int k) {
  if(len1 != len2 || kord1 == NULL || kord2 == NULL) { return(-1.0); } // Exit if different lengths
  size_t klen = len1 - k + 1;
  uint16_t dotsum = 0;
  for(size_t i=0; i<klen; i++) {
    dotsum += (kord1[i] == kord2[i]);
  }
  double dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}

#if DADAS_SIMD

/* ---------- Mínimos de vectores de conteo de k-mers, 16 bits ---------- */
// 8 carriles de 16 bits, suma con desborde (igual que _mm_add_epi16); los conteos son < 2^15,
// así que el mínimo con o sin signo es el mismo.
static inline double kmer_dist_simd16(const uint16_t *kv1, int len1, const uint16_t *kv2, int len2, int k) {
  size_t n_kmer = 1 << (2*k); // asume n_kmer divisible por 8
  uint16_t dst[8];
#if defined(__x86_64__) || defined(__x86_64)
  __m128i vsum = _mm_set1_epi16(0);
  for(const uint16_t *end = kv1 + n_kmer; kv1 < end; kv1 += 8, kv2 += 8) {
    __m128i kk1 = _mm_loadu_si128((const __m128i*) kv1);
    __m128i kk2 = _mm_loadu_si128((const __m128i*) kv2);
    vsum = _mm_add_epi16(vsum, _mm_min_epi16(kk1, kk2));
  }
  _mm_storeu_si128((__m128i*) dst, vsum);
#else
  uint16x8_t vsum = vdupq_n_u16(0);
  for(const uint16_t *end = kv1 + n_kmer; kv1 < end; kv1 += 8, kv2 += 8) {
    vsum = vaddq_u16(vsum, vminq_u16(vld1q_u16(kv1), vld1q_u16(kv2)));
  }
  vst1q_u16(dst, vsum);
#endif
  uint16_t dotsum = 0;
  for(int i=0; i<8; i++) { dotsum += dst[i]; }
  double dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}

/* ---------- Mínimos de vectores de conteo de k-mers, 8 bits ---------- */
// 16 carriles de 8 bits con suma saturada (igual que _mm_adds_epu8). Si algún carril satura (255),
// devuelve -1 y quien llama repite con la versión de 16 bits.
static inline double kmer_dist_simd8(const uint8_t *kv1, int len1, const uint8_t *kv2, int len2, int k) {
  size_t n_kmer = 1 << (2*k); // asume n_kmer divisible por 16
  uint8_t dst[16];
#if defined(__x86_64__) || defined(__x86_64)
  __m128i vsum = _mm_set1_epi8(0x00);
  for(const uint8_t *end = kv1 + n_kmer; kv1 < end; kv1 += 16, kv2 += 16) {
    __m128i kk1 = _mm_loadu_si128((const __m128i*) kv1);
    __m128i kk2 = _mm_loadu_si128((const __m128i*) kv2);
    vsum = _mm_adds_epu8(vsum, _mm_min_epu8(kk1, kk2));
  }
  _mm_storeu_si128((__m128i*) dst, vsum);
#else
  uint8x16_t vsum = vdupq_n_u8(0);
  for(const uint8_t *end = kv1 + n_kmer; kv1 < end; kv1 += 16, kv2 += 16) {
    vsum = vqaddq_u8(vsum, vminq_u8(vld1q_u8(kv1), vld1q_u8(kv2)));
  }
  vst1q_u8(dst, vsum);
#endif
  uint16_t dotsum = 0;
  bool overflow = false;
  for(int i=0; i<16; i++) {
    if(dst[i] == 255) { overflow = true; }
    dotsum += dst[i];
  }
  if(overflow) { return(-1.); }
  double dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}

/* ---------- k-mers en orden (para detectar alineamientos sin gaps) ---------- */
static inline double kord_dist_simd(const uint16_t *kord1, int len1, const uint16_t *kord2, int len2, int k) {
  if(kord1 == NULL || kord2 == NULL) { return(-1.0); }
  size_t klen = (len1 < len2 ? len1 : len2) - k + 1;
  size_t n_vec = klen - (klen % 8); // vectoriza en bloques de 8
  uint16_t dst[8];
  const uint16_t *end = kord1 + n_vec;
#if defined(__x86_64__) || defined(__x86_64)
  __m128i vsum = _mm_set1_epi16(0);
  for(; kord1 < end; kord1 += 8, kord2 += 8) {
    __m128i keq = _mm_cmpeq_epi16(_mm_loadu_si128((const __m128i*) kord1), _mm_loadu_si128((const __m128i*) kord2));
    vsum = _mm_sub_epi16(vsum, keq); // 0xFFFF = -1
  }
  _mm_storeu_si128((__m128i*) dst, vsum);
#else
  uint16x8_t vsum = vdupq_n_u16(0);
  for(; kord1 < end; kord1 += 8, kord2 += 8) {
    vsum = vsubq_u16(vsum, vceqq_u16(vld1q_u16(kord1), vld1q_u16(kord2))); // 0xFFFF = -1
  }
  vst1q_u16(dst, vsum);
#endif
  uint16_t dotsum = 0;
  for(int i=0; i<8; i++) { dotsum += dst[i]; }
  for(size_t i=n_vec; i<klen; i++, kord1++, kord2++) { // resto, escalar
    if(*kord1 == *kord2) { dotsum++; }
  }
  double dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}

#endif // DADAS_SIMD
#endif // DADAS_KMERS_SIMD_H
