#include <Rcpp.h>
#include "dada.h"
using namespace Rcpp;

#define MIN(a,b) (((a)<(b))?(a):(b))

// DadaS: R builds packages with -O2. With GCC < 14 that does NOT vectorize these loops (the "very-cheap"
// cost model of -O2 rejects loops that need an epilogue), so the "vectorized" aligner ran scalar on Linux.
// We force vectorization for this file only. Integer arithmetic: identical results.
// (clang, the macOS compiler, already vectorizes at -O2; there this changes nothing.)
#if defined(__GNUC__) && !defined(__clang__)
#define DADAS_VECTORIZE __attribute__((optimize("tree-vectorize", "vect-cost-model=dynamic")))
#else
#define DADAS_VECTORIZE
#endif

// DadaS: one anti-diagonal loop with the diagonal-score computation fused in.
// Previously diag_buf was filled in a separate scalar loop (s1 walked backwards, s2 forwards: not
// vectorizable) and read here. Now s1r is s1 reversed, so a[] and b[] advance together and the whole
// loop vectorizes. Same integer operations and same tie-breaking: identical results.
//   SWAP=false: precedence up > left > diag      SWAP=true: left > up > diag
template<bool SWAP>
DADAS_VECTORIZE
static inline void dploop_fused(const int16_t *__restrict__ ptr_left, const int16_t *__restrict__ ptr_up,
                                const int16_t *__restrict__ ptr_d2, const char *__restrict__ a, const char *__restrict__ b,
                                int16_t *__restrict__ d, int16_t *__restrict__ p,
                                int16_t match, int16_t mismatch, int16_t gap_p, size_t n) {
  for(size_t k=0; k<n; k++) {
    int16_t left = ptr_left[k] + gap_p;
    int16_t up = ptr_up[k] + gap_p;
    int16_t diag = ptr_d2[k] + (a[k] == b[k] ? match : mismatch);
    int16_t entry, pentry;
    if(SWAP) {
      entry = left >= up ? left : up;
      pentry = left >= up ? 2 : 3;
    } else {
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
    }
    pentry = entry >= diag ? pentry : 1;
    entry = entry >= diag ? entry : diag;
    d[k] = entry;
    p[k] = pentry;
  }
}

void parr(int16_t *arr, int nrow, int ncol) {
  int col, row;
  for(row=0;row<nrow;row++) {
    for(col=0;col<ncol;col++) {
      Rprintf("%05d\t", arr[row*ncol+col]);
    }
    Rprintf("\n");
  }
}


static void nw_vec_fill(const char *&s1_io, size_t &len1_io, const char *&s2_io, size_t &len2_io, int16_t match, int16_t mismatch, int16_t gap_p, int16_t end_gap_p, int band, NWWork &w) {
  // local copies (written back at the end): references would force reloads inside the hot loops
  const char *s1 = s1_io, *s2 = s2_io;
  size_t len1 = len1_io, len2 = len2_io;
  size_t row, col, ncol, nrow, foo;
  size_t i,j;
  int16_t d_free;
  size_t start_col;
  size_t col_min, col_max, even;
  size_t i_max, j_min;
  int16_t *ptr_left, *ptr_up, *ptr_d, *ptr_p;
  bool swap = false;
  bool recalc_left = false, recalc_right = false;
  const char *ptr_const_char;
  char *ptr_char;

  if(len1 > len2) { // Ensure s1 is the shorter sequences
    ptr_const_char = s1;
    s1 = s2;
    s2 = ptr_const_char;
    swap = true;
    foo = len1;
    len1 = len2;
    len2 = foo;
  }
  if(band < 0) { band = len2; }
  
  // Allocate the DP matrices
  start_col = 1 + (1+(band<len1 ? band : len1))/2;
  //  ncol = 3 + (len1+len2+1)/2; // 3 = left boundary + center + right boundary !!!
  ncol = 2 + start_col + ((len2-len1+band)<len2 ? (len2-len1+band) : len2)/2;
  nrow = len1 + len2 + 1;
  // DadaS: per-thread buffers reused across calls instead of a malloc/free of each DP matrix per alignment.
  // Every cell that is read is written first (as before, the matrices were not initialized either).
  if(w.d.size() < ncol * nrow) { w.d.resize(ncol * nrow); w.p.resize(ncol * nrow); }
  if(w.s1r.size() < len1 + 1) { w.s1r.resize(len1 + 1); }
  int16_t *d = w.d.data();
  int16_t *p = w.p.data();
  char *s1r = w.s1r.data(); // s1 reversed (see dploop_fused)
  for(i=0;i<len1;i++) { s1r[i] = s1[len1-1-i]; }
  
  // For banding issues later on
  int16_t fill_val = INT16_MIN - MIN(MIN(mismatch, gap_p), MIN(match, 0));
  for(row=0;row<nrow;row++) {
    d[row*ncol] = fill_val;
    d[row*ncol+1] = fill_val;
    d[row*ncol+ncol-2] = fill_val;
    d[row*ncol+ncol-1] = fill_val;
  }

  // Fill out starting point
  d[start_col] = 0;
  p[start_col] = 0; // Should never be queried
  
  // Fill out "left" "column" of d, p.
  row=1;
  col=start_col-1;
  ptr_d=&d[ncol]; // start of 1st row
  ptr_p=&p[ncol];
  d_free = end_gap_p; // end-gap score
  while(row < (1 + (band < len1 ? band : len1))) {
    ptr_d[col] = d_free; // end gap
    ptr_p[col] = 3;
    if(row%2==0) {
      col--;
    }
    row++;
    d_free += end_gap_p;
    ptr_d += ncol;
    ptr_p += ncol;
  }
  
  // Fill out "top" "row" of d, p.
  row=1;
  col=start_col;
  ptr_d=&d[ncol]; // start of 1st row
  ptr_p=&p[ncol];
  d_free = end_gap_p; // end-gap score
  while(row < (1+(band+len2-len1 < len2 ? band+len2-len1 : len2))) {
    ptr_d[col] = d_free;
    ptr_p[col] = 2;
    if(row%2 == 1) {
      col++;
    }
    row++;
    d_free += end_gap_p;
    ptr_d += ncol;
    ptr_p += ncol;
  }

  // Fill out DP matrix (Row 0/1 taken care of by ends-free)
  row = 2;
  col_min = start_col; // Do not fill out the ends-free cells
  col_max = start_col;
  i_max = 0; // 1st nt
  j_min = 0; // 1st nt
  even = 1; // TRUE
  
  while(row <= (len1+len2)) {
      // Fill out row
    // s1[i_max - k] == s1r[len1-1-i_max + k]  and  s2[j_min + k], for k = 0..n-1
    ptr_left = &d[(row-1)*ncol + col_min-even];
    ptr_up = &d[(row-1)*ncol + col_min+1-even];
    if(swap) {
      dploop_fused<true>(ptr_left, ptr_up, &d[(row-2)*ncol + col_min], &s1r[len1-1-i_max], &s2[j_min],
                         &d[row*ncol + col_min], &p[row*ncol + col_min], match, mismatch, gap_p, col_max-col_min+1);
    } else {
      dploop_fused<false>(ptr_left, ptr_up, &d[(row-2)*ncol + col_min], &s1r[len1-1-i_max], &s2[j_min],
                          &d[row*ncol + col_min], &p[row*ncol + col_min], match, mismatch, gap_p, col_max-col_min+1);
    }
    
    // Offset to the boundary after top wedge (w/ prefilled boundary) is done
    if(row==(band<len1 ? band : len1)) { 
      col_min--;
      i_max++;
      j_min--;
    }
    if(row == (band+len2-len1 < len2 ? band+len2-len1 : len2)) {
      col_max++;
    }
    
    // Recalculate ends-free cells for lower boundary
    if(end_gap_p > gap_p) {
      // Left column
      if(recalc_left) { // past first row of lower tri
        d_free = ptr_left[0] + end_gap_p; // first cell of left array
        if(d_free > d[row*ncol + col_min]) { // ends-free gap is better
          d[row*ncol + col_min] = d_free;
          p[row*ncol + col_min] = 2;
        } else if(!swap && d_free == d[row*ncol + col_min] && p[row*ncol + col_min] == 1) { // left gap takes precedence over diagonal move (for consistency)
          p[row*ncol + col_min] = 2;
        } else if(swap && d_free == d[row*ncol + col_min] && p[row*ncol + col_min] != 2) { // left-is-up, and takes precedence other moves
          p[row*ncol + col_min] = 2;
        }
      }
      if(i_max == len1-1) { recalc_left = true; }
      // Right column
      if(recalc_right) {
        d_free = ptr_up[col_max-col_min] + end_gap_p; // last cell of up array
        if(d_free > d[row*ncol + col_max]) { // ends-free gap is better
          d[row*ncol + col_max] = d_free;
          p[row*ncol + col_max] = 3;
        } else if(!swap && d_free == d[row*ncol + col_max] && p[row*ncol + col_max] != 3) { // up gap takes precedence over left or diagonal move (for consistency)
          p[row*ncol + col_max] = 3;
        } else if(swap && d_free == d[row*ncol + col_max] && p[row*ncol + col_max] == 1) { // up-is-left, and takes precedence over diagonal
          p[row*ncol + col_max] = 3;
        }
      }
      if((row+1)/2 + col_max - start_col == len2) { recalc_right = true; }
      // j_max (1-index) = (row+1)/2 + col_max - start_col
    }
    // Update the indices
    if(row < band && row < len1) { // upper tri for seq1
      if(even) { col_min--; }
      i_max++;
    } else if(i_max < len1-1) { // banded area
      if(band%2 == 0) {
        if(even) { j_min++; }
        else { i_max++; }
      } else { // odd band
        if(even) { col_min--; i_max++; }
        else { col_min++; j_min++; }
      }
    } else { // lower tri for seq1
      if(!even) { col_min++; }
      j_min++;
    }

    if(row<(band+len2-len1 < len2 ? band+len2-len1 : len2)) {
      if(!even) { col_max++; }
    } else if((row+1)/2 + col_max - start_col < len2) { // "j_max" (1-index) < len2
      if((band+len2-len1) % 2 == 0) { // even band (including the extra band from length difference)
        if(even) { col_max--; }
        else { col_max++; }
      } // no action on odd band
    } else {
      if(even) { col_max--; }
    }
    
    row++;
    even = 1 - even;
  }

  w.ncol = ncol;
  w.start_col = start_col;
  w.swap = swap;
  s1_io = s1; s2_io = s2; len1_io = len1; len2_io = len2;
}

char **nwalign_vectorized2(const char *s1, size_t len1, const char *s2, size_t len2, int16_t match, int16_t mismatch, int16_t gap_p, int16_t end_gap_p, int band) {
  size_t i, j;
  char *ptr_char;
  NWWork w; // local buffers (this path is used outside the denoising loop)
  nw_vec_fill(s1, len1, s2, len2, match, mismatch, gap_p, end_gap_p, band, w);
  const int16_t *p = w.p.data();
  size_t ncol = w.ncol, start_col = w.start_col, nrow = len1 + len2 + 1;
  bool swap = w.swap;

//  if(debug) {
//    Rprintf("D MATRIX:.....\n");
//    parr(d, len1+len2+1, ncol);
//    Rprintf("\n");
//    Rprintf("Score: %d\n", d[(len1+len2)*ncol + (2*start_col+len2-len1)/2]);
//  }

  char *al0 = (char *) malloc((nrow+1) * sizeof(char));
  char *al1 = (char *) malloc((nrow+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) dada_stop("Memory allocation failed.");
  
  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;
  
  while ( i > 0 || j > 0 ) {
    switch ( p[(i+j)*ncol + (2*start_col+j-i)/2] ) {
      case 1:
        al0[len_al] = s1[--i];
        al1[len_al] = s2[--j];
        break;
      case 2:
        al0[len_al] = '-';
        al1[len_al] = s2[--j];
        break;
      case 3:
        al0[len_al] = s1[--i];
        al1[len_al] = '-';
        break;
      default:
        dada_stop("N-W Align out of range: len1/2=(%i, %i), nrow,ncol=(%i,%i), ij=(%i,%i).", (int) len1, (int) len2, (int) nrow, (int) ncol, (int) i, (int) j);
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';
  
  
  // Return to input ordering
  if(swap) {
    ptr_char = al0;
    al0 = al1;
    al1 = ptr_char;
  }
  
  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  dada_stop("Failed memory allocation.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  dada_stop("Failed memory allocation.");
  
  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
  free(al0);
  free(al1);
  
  return al;
}


// DadaS: same alignment as nwalign_vectorized2() followed by al2subs(), without building the two aligned
// strings. The traceback writes the center->raw position map, the substitutions and the internal indel count
// directly. s1 must be the center (reference) and s2 the raw, as in sub_new(). Results are identical to
// al2subs(nwalign_vectorized2(...)): checked pair by pair with C_check_fused_subs() on the test data.
Sub *nwalign_vectorized2_sub(const char *s1, size_t len1, const char *s2, size_t len2, int16_t match, int16_t mismatch, int16_t gap_p, int16_t end_gap_p, int band, NWWork *wp) {
  const char *c = s1, *r = s2;         // center and raw, before any swap
  size_t lenc = len1;
  NWWork wlocal;
  NWWork &w = wp ? *wp : wlocal;       // caller's reusable buffers, or temporary ones
  nw_vec_fill(s1, len1, s2, len2, match, mismatch, gap_p, end_gap_p, band, w);
  const int16_t *p = w.p.data();
  size_t ncol = w.ncol, start_col = w.start_col;
  bool swap = w.swap;

  // Traceback (same moves as nwalign_vectorized2). Columns are recorded from the end of the alignment.
  size_t maxcols = len1 + len2;
  if(w.colc.size() < maxcols) { w.colc.resize(maxcols); w.colr.resize(maxcols); }
  int *cc = w.colc.data(), *cr = w.colr.data(); // center / raw position of each column, -1 = gap
  size_t ncolumns = 0;
  long i = len1, j = len2, a, b;
  while ( i > 0 || j > 0 ) {
    switch ( p[(i+j)*ncol + (2*start_col+j-i)/2] ) {
      case 1: a = --i; b = --j; break;
      case 2: a = -1;  b = --j; break;
      case 3: a = --i; b = -1;  break;
      default:
        dada_stop("N-W Align out of range: len1/2=(%i, %i), ij=(%i,%i).", (int) len1, (int) len2, (int) i, (int) j);
    }
    if(swap) { cc[ncolumns] = (int) b; cr[ncolumns] = (int) a; } // s1/s2 were swapped: s1 is the raw
    else     { cc[ncolumns] = (int) a; cr[ncolumns] = (int) b; }
    ncolumns++;
  }

  // Count substitutions; find the first and last columns where both sequences have a nucleotide
  unsigned int nsubs = 0;
  long first = -1, last = -1; // in forward column order
  for(size_t k=0; k<ncolumns; k++) {
    size_t f = ncolumns - 1 - k; // forward index of this recorded column
    int pc = cc[k], pr = cr[k];
    if(pc >= 0 && pr >= 0) {
      if(c[pc] != r[pr] && c[pc] != 5 && r[pr] != 5) nsubs++;
      if(last < 0) last = (long) f;
      first = (long) f;
    }
  }
  Sub *sub = (Sub *) malloc(sizeof(Sub)); //E
  if (sub == NULL)  dada_stop("Memory allocation failed.");
  sub->len0 = lenc;
  sub->map = (uint16_t *) malloc(lenc * sizeof(uint16_t)); //E
  sub->pos = (uint16_t *) malloc(nsubs * sizeof(uint16_t)); //E
  sub->nt0 = (char *) malloc(nsubs); //E
  sub->nt1 = (char *) malloc(nsubs); //E
  if (sub->map == NULL || sub->pos == NULL || sub->nt0 == NULL || sub->nt1 == NULL) {
    dada_stop("Memory allocation failed.");
  }
  // Forward pass: map, substitutions in center order, internal indel events (gap runs between first and last)
  unsigned int s = 0, nindel = 0;
  bool prev_gap = false;
  for(size_t f=0; f<ncolumns; f++) {
    size_t k = ncolumns - 1 - f;
    int pc = cc[k], pr = cr[k];
    if(pc >= 0) {
      sub->map[pc] = (pr >= 0) ? (uint16_t) pr : GAP_GLYPH;
      if(pr >= 0 && c[pc] != r[pr] && c[pc] != 5 && r[pr] != 5) {
        sub->pos[s] = (uint16_t) pc; sub->nt0[s] = c[pc]; sub->nt1[s] = r[pr]; s++;
      }
    }
    if(first >= 0 && (long) f >= first && (long) f <= last) {
      bool gap = (pc < 0) || (pr < 0);
      if(gap && !prev_gap) { nindel++; }
      prev_gap = gap;
    }
  }
  sub->nsubs = nsubs;
  sub->nindel = nindel;
  sub->q0 = NULL;
  sub->q1 = NULL;
  return sub;
}

// [[Rcpp::export]]
Rcpp::CharacterVector C_nwvec(std::vector<std::string> s1, std::vector<std::string> s2, int16_t match, int16_t mismatch, int16_t gap_p, int band, bool endsfree) {
  char **al;
  int i;
  if(s1.size() != s2.size()) {
    dada_stop("Character vectors to be aligned must be of equal length.");
  }
  Rcpp::CharacterVector rval(s1.size()*2);
  
  for(i=0;i<s1.size();i++) {
    if(endsfree) {
      al = nwalign_vectorized2(s1[i].c_str(), s1[i].length(), s2[i].c_str(), s2[i].length(), match, mismatch, gap_p, 0, (size_t) band);
    } else {
      al = nwalign_vectorized2(s1[i].c_str(), s1[i].length(), s2[i].c_str(), s2[i].length(), match, mismatch, gap_p, gap_p, (size_t) band);
    }

    rval[2*i] = std::string(al[0]);
    rval[2*i+1] = std::string(al[1]);
    free(al[0]);
    free(al[1]);
    free(al);
  }
  return(rval);
}
