#include "dada.h"
#include "kmers_simd.h"
#include <Rcpp.h>
///#ifdef _WIN32 //  Windows
/*  //do arm stuff
#endif
#if defined(_WIN32) && !defined(__MINGW32__) //  Windows and not MINGW
#define cpuid(info, x)    __cpuidex(info, x, 0)
#else //  GCC Intrinsics
#include <cpuid.h>
void cpuid(int info[4], int InfoType){
  __cpuid_count(InfoType, 0, info[0], info[1], info[2], info[3]);
}
#endif */

using namespace Rcpp;
//' @useDynLib dadaS
//' @importFrom Rcpp evalCpp

void run_dada(B *bb, Rcpp::NumericMatrix errMat, 
              int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, 
              bool detect_singletons, int max_clust, double min_fold, int min_hamming, int min_abund, 
              bool vectorized_alignment, bool multithread, bool verbose, 
              int SSE, bool gapless, bool greedy, AlignCache *cache, const double *emax);

//------------------------------------------------------------------
// C interface to run DADA on the provided unique sequences/abundance pairs. 
// 
// [[Rcpp::export]]
Rcpp::List dada_uniques(std::vector< std::string > seqs, std::vector<int> abundances, std::vector<bool> priors,
                        Rcpp::NumericMatrix err,
                        Rcpp::NumericMatrix quals,
                        int match, int mismatch, int gap,
                        bool use_kmers, double kdist_cutoff,
                        int band_size,
                        double omegaA, double omegaP, double omegaC, bool detect_singletons,
                        int max_clust,
                        double min_fold, int min_hamming, int min_abund,
                        bool use_quals,
                        bool final_consensus,
                        bool vectorized_alignment,
                        int homo_gap,
                        bool multithread,
                        bool verbose,
                        int SSE,
                        bool gapless,
                        bool greedy,
                        double indel_rate,
                        SEXP align_cache,
                        std::vector<int> qh_values,
                        std::vector<int> qh_index,
                        std::vector<int> qh_counts) {

  unsigned int i, r, index, pos, nraw, maxlen, minlen;
  Raw *raw;
  
  /********** INPUT VALIDATION *********/
  // Check lengths of seqs and abundances vectors
  nraw = seqs.size();
  if(nraw == 0) { Rcpp::stop("Zero input sequences."); }
  if(nraw != abundances.size()) { Rcpp::stop("Sequence and abundance vectors had different lengths."); }
  if(nraw != priors.size()) { Rcpp::stop("Sequence and priors vectors had different lengths."); }
  maxlen=0;
  minlen=SEQLEN;
  for(index=0;index<nraw;index++) {
    if(seqs[index].length() > maxlen) { maxlen = seqs[index].length(); }
    if(seqs[index].length() < minlen) { minlen = seqs[index].length(); }
  }
  if(maxlen >= SEQLEN) { Rcpp::stop("Input sequences exceed the maximum allowed string length."); }
  if(minlen <= KMER_SIZE) { Rcpp::stop("Input sequences must all be longer than the kmer-size (%i).", KMER_SIZE); }
  
  // Check for presence of quality scores and their lengths
  bool has_quals = false;
  if(quals.nrow() > 0) { // Each sequence is a COLUMN, each row is a POSITION
    has_quals = true;
    if(quals.nrow() != maxlen) {
      Rcpp::stop("Sequence must have associated qualities for each nucleotide position.");
    }
  }
  // Check error matrix
  if(indel_rate >= 1) { Rcpp::stop("INDEL_RATE must be < 1."); }
  if(err.nrow() != 16) {
    Rcpp::stop("Error matrix must have 16 rows.");
  }
  // Check for SSE2+ support
  // Code adapted from Stack Overflow (Mysticial) 
  // https://stackoverflow.com/questions/6121792/how-to-check-if-a-cpu-supports-the-sse3-instruction-set
/*  bool HW_SSE = false;
  bool HW_SSE2 = false;
  bool HW_SSE3 = false;
  int info[4];
  cpuid(info, 0);
  int nIds = info[0];

  //  Detect Features
  if (nIds >= 0x00000001){
    cpuid(info,0x00000001);
    HW_SSE    = (info[3] & ((int)1 << 25)) != 0;
    HW_SSE2   = (info[3] & ((int)1 << 26)) != 0;
    HW_SSE3   = (info[2] & ((int)1 <<  0)) != 0;
  }
  if(!(HW_SSE && HW_SSE2)) { 
    Rprintf("SSE2 not supported.\n");
    SSE = 0;
  } */
  if(!DADAS_SIMD) { SSE=0; } // DadaS: SSE2 on x86 or NEON on ARM (see kmers_simd.h)
  
  /********** CONSTRUCT RAWS *********/
  // DadaS: every allocation below is released in cleanup(), which also runs when an exception leaves this
  // function (an error in a worker thread, or Ctrl-C via Rcpp::checkUserInterrupt). Previously an interrupted
  // or failed dada() left all of this memory allocated for the rest of the R session.
  char seq[SEQLEN];
  double qual[SEQLEN];
  Raw **raws = NULL;
  uint8_t *k8 = NULL;
  uint16_t *k16 = NULL;
  uint16_t *kord = NULL;
  double *emax = NULL;
  B *bb = NULL;
  Sub **subs = NULL;
  Sub **birth_subs = NULL;
  unsigned int nclust_subs = 0;
  auto cleanup = [&]() {
    if(subs) { for(index=0;index<nraw;index++) { sub_free(subs[index]); } free(subs); subs = NULL; }
    if(birth_subs) { for(i=0;i<nclust_subs;i++) { sub_free(birth_subs[i]); } free(birth_subs); birth_subs = NULL; }
    if(bb) { b_free(bb); bb = NULL; }
    if(raws) { for(index=0;index<nraw;index++) { if(raws[index]) raw_free(raws[index]); } free(raws); raws = NULL; }
    free(k8); k8 = NULL; free(k16); k16 = NULL; free(kord); kord = NULL; free(emax); emax = NULL;
  };

  // DadaS: optional cache of alignment results shared by successive dada_uniques calls on the same sample
  // (the self-consistency rounds of learnErrors). See AlignCache in dada.h.
  AlignCache *cache = NULL;
  if(!Rf_isNull(align_cache)) {
    Rcpp::XPtr<AlignCache> xp(align_cache);
    cache = xp.get();
    if(cache && cache->nraw != nraw) { cache->by_center.clear(); cache->nraw = nraw; }
  }

  Rcpp::List rval;
  try {
    raws = (Raw **) calloc(nraw, sizeof(Raw *)); //E
    if (raws == NULL)  dada_stop("Memory allocation failed.");
    // Construct a raw for each input sequence, store in raws[index]
    for (index = 0; index < nraw; index++) {
      strcpy(seq, seqs[index].c_str());
      nt2int(seq, seq);
      if(has_quals) {
        for(pos=0;pos<seqs[index].length();pos++) {
          qual[pos] = quals(pos, index);
        }
        raws[index] = raw_new(seq, qual, abundances[index], priors[index]);
      } else {
        raws[index] = raw_new(seq, NULL, abundances[index], priors[index]);
      }
      raws[index]->index = index;
    }
    
    // DadaS (PER_READ_QUALS): attach the per-read quality histogram to the uniques that have one
    if(qh_values.size() > 0) {
      if(!use_quals || !has_quals) dada_stop("PER_READ_QUALS requires quality scores.");
      unsigned int V = qh_values.size(), nh = qh_index.size();
      if(qh_counts.size() != (size_t) nh * maxlen * V) dada_stop("Quality histogram has the wrong size.");
      for(unsigned int v=0; v<V; v++) { if(qh_values[v] < 0 || qh_values[v] >= err.ncol()) dada_stop("Quality score %i outside the error matrix.", qh_values[v]); }
      for(unsigned int k=0; k<nh; k++) {
        if(qh_index[k] < 0 || (unsigned int) qh_index[k] >= nraw) dada_stop("Quality histogram index out of range.");
        raw = raws[qh_index[k]];
        raw->nqv = V;
        raw->qvals = qh_values.data();
        raw->qcount = (uint32_t *) malloc(raw->length * V * sizeof(uint32_t));
        raw->qfrac = (double *) malloc(raw->length * V * sizeof(double));
        if(!raw->qcount || !raw->qfrac) dada_stop("Memory allocation failed.");
        for(pos=0; pos<raw->length; pos++) {
          unsigned int tot = 0;
          for(unsigned int v=0; v<V; v++) {
            int c = qh_counts[k + (size_t) nh * (pos + (size_t) maxlen * v)];
            raw->qcount[pos*V + v] = (uint32_t) c;
            raw->qfrac[pos*V + v] = ((double) c) / raw->reads;
            tot += c;
          }
          if(tot != raw->reads) dada_stop("Quality histogram does not match the abundance of unique %u.", qh_index[k] + 1);
        }
      }
    }

    if(use_kmers) {
      /// Add uint8_t kmer index in contiguous memory block
      size_t n_kmer = 1 << (2*KMER_SIZE);
      k8 = (uint8_t *) malloc(nraw * n_kmer * sizeof(uint8_t)); //E
      if (k8 == NULL)  dada_stop("Memory allocation failed.");
      for(index=0;index<nraw;index++) {
        raw = raws[index];
        raw->kmer8 = &k8[index*n_kmer];
        assign_kmer8(raw->kmer8, raw->seq, KMER_SIZE);
      }
      /// Add uint16_t kmer index in contiguous memory block
      k16 = (uint16_t *) malloc(nraw * n_kmer * sizeof(uint16_t)); //E
      if (k16 == NULL)  dada_stop("Memory allocation failed.");
      for(index=0;index<nraw;index++) {
        raw = raws[index];
        raw->kmer = &k16[index*n_kmer];
        assign_kmer(raw->kmer, raw->seq, KMER_SIZE);
      }
      /// Add uint16_t ordered kmer record in contiguous memory block
      kord = (uint16_t *) malloc(nraw * maxlen * sizeof(uint16_t)); //E
      if (kord == NULL)  dada_stop("Memory allocation failed.");
      for(index=0;index<nraw;index++) {
        raw = raws[index];
        raw->kord = &kord[index*maxlen];
        assign_kmer_order(raw->kord, raw->seq, KMER_SIZE);
      }
    } else {
      for(index=0;index<nraw;index++) {
        raw = raws[index];
        raw->kmer8 = NULL;
        raw->kmer = NULL;
        raw->kord = NULL;
      }
    }

    // DadaS: per-unique upper bound on the error rate of any substitution it could carry (the largest
    // off-diagonal entry of the error matrix at any of its quality scores). Used with the alignment cache.
    if(cache) {
      unsigned int ncol = err.ncol();
      emax = (double *) malloc(nraw * sizeof(double));
      if(emax == NULL) dada_stop("Memory allocation failed.");
      std::vector<double> colmax(ncol, 0.0);
      for(unsigned int col=0; col<ncol; col++) {
        for(unsigned int t=0; t<16; t++) {
          if(t % 5 == 0) continue; // A2A, C2C, G2G, T2T
          if(err(t, col) > colmax[col]) colmax[col] = err(t, col);
        }
      }
      for(index=0;index<nraw;index++) {
        raw = raws[index];
        double m = 0.0;
        if(use_quals && raw->qcount) { // DadaS: any quality score present at a position can carry a substitution
          for(pos=0;pos<raw->length && m < 1.0;pos++) {
            for(unsigned int v=0; v<raw->nqv; v++) {
              if(raw->qcount[pos*raw->nqv + v] == 0) continue;
              unsigned int q = raw->qvals[v];
              if(q >= ncol) { m = 1.0; break; }
              if(colmax[q] > m) m = colmax[q];
            }
          }
        } else if(use_quals && raw->qual) {
          for(pos=0;pos<raw->length;pos++) {
            unsigned int q = raw->qual[pos];
            if(q >= ncol) { m = 1.0; break; } // out of range: compute_lambda will stop; no pruning
            if(colmax[q] > m) m = colmax[q];
          }
        } else {
          m = colmax[0];
        }
        emax[index] = m > 1.0 ? 1.0 : m;
      }
    }

    /********** RUN DADA *********/
    bb = b_new(raws, nraw, omegaA, omegaP, use_quals); // New cluster with all sequences in 1 bi
    bb->indel_rate = indel_rate; // DadaS: <= 0 reproduces DADA2 exactly (indels ignored)
    run_dada(bb, err, match, mismatch, gap, homo_gap, use_kmers, kdist_cutoff, band_size, 
             detect_singletons, max_clust, min_fold, min_hamming, min_abund, 
             vectorized_alignment, multithread, verbose, SSE, gapless, greedy, cache, emax);

    /********** MAKE OUTPUT *********/
    // Create subs for all the relevant alignments
    subs = (Sub **) calloc(bb->nraw, sizeof(Sub *)); //E
    nclust_subs = bb->nclust;
    birth_subs = (Sub **) calloc(bb->nclust, sizeof(Sub *)); //E
    if(!subs || !birth_subs) dada_stop("Memory allocation failed.");
    
    // Parallelize over clusters
    struct FinalSubsParallel : public RcppParallel::Worker
    {
      // source data
      B *b;
      // destination Sub arrays
      Sub **subs;
      Sub **birth_subs;
      // parameters
      int match, mismatch, gap, homo_gap;
      int band_size;
      bool use_kmers, vectorized_alignment;
      int SSE;
      bool gapless;
      
      // initialize with source and destination
      FinalSubsParallel(B *b, Sub **subs, Sub **birth_subs,
                      int match, int mismatch, int gap, int homo_gap,
                      int band_size, bool use_kmers, bool vectorized_alignment, int SSE, bool gapless) 
        : b(b), subs(subs), birth_subs(birth_subs), match(match), mismatch(mismatch), gap(gap), homo_gap(homo_gap), 
          band_size(band_size), use_kmers(use_kmers), vectorized_alignment(vectorized_alignment), SSE(SSE), gapless(gapless) {}
      
      // Perform sequence comparison
      void operator()(std::size_t begin, std::size_t end) {
        Raw *raw;
        NWWork nww; // DadaS: reused aligner buffers
        for(std::size_t i=begin;i<end;i++) {
          for(unsigned int r=0;r<b->bi[i]->nraw;r++) {
            raw = b->bi[i]->raw[r];
            subs[raw->index] = sub_new(b->bi[i]->center, raw, match, mismatch, gap, homo_gap, false, 1.0, band_size, vectorized_alignment, SSE, gapless, &nww);
          }
          // Make birth sub for that cluster
          if(i==0) { birth_subs[i] = NULL; }
          else {
            birth_subs[i] = sub_new(b->bi[b->bi[i]->birth_comp.i]->center, b->bi[i]->center, match, mismatch, gap, homo_gap, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless);
          }
        } // for(std::size_t i=begin;i<end;i++)
      }
    };
    
    if(multithread) {
      FinalSubsParallel finalSubsParallel(bb, subs, birth_subs, match, mismatch, gap, homo_gap, band_size, use_kmers, vectorized_alignment, SSE, gapless);
      RcppParallel::parallelFor(0, bb->nclust, finalSubsParallel, GRAIN_SIZE);
    } else { // Non-Parallel implementation
      NWWork nww; // DadaS: reused aligner buffers
      for(i=0;i<bb->nclust;i++) {
        // Make subs for members of that cluster
        for(r=0;r<bb->bi[i]->nraw;r++) {
          raw = bb->bi[i]->raw[r];
          subs[raw->index] = sub_new(bb->bi[i]->center, raw, match, mismatch, gap, homo_gap, false, 1.0, band_size, vectorized_alignment, SSE, gapless, &nww);
        }
        // Make birth sub for that cluster
        if(i==0) { birth_subs[i] = NULL; }
        else {
          birth_subs[i] = sub_new(bb->bi[bb->bi[i]->birth_comp.i]->center, bb->bi[i]->center, match, mismatch, gap, homo_gap, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless);
        }
      }
    }

    // Assign final w/in cluster P
    Rcpp::NumericVector Rpraw(nraw); ///! temporary output for evaluation purposes
    for(i=0;i<bb->nclust;i++) {
      for(r=0;r<bb->bi[i]->nraw;r++) {
        raw = bb->bi[i]->raw[r];
        index = raw->index;
        if(bb->bi[i]->center == raw) { // centers always have correct=true
          raw->p = 1.0; // Due to fp precision, possible to get 1's from non-centers in below calc
        } else {
          raw->p = calc_pA(raw->reads, raw->comp.lambda * bb->bi[i]->reads, true);
          if(raw->p < omegaC) { raw->correct = false; } // Don't correct!
        }
        Rpraw[index] = raw->p; ///! temporary output for evaluation purposes
      }
    }
    
    Rcpp::DataFrame df_clustering = b_make_clustering_df(bb, subs, birth_subs, has_quals);
    // DadaS: with USE_QUALS=FALSE the error model has one column, so transitions are not split by quality.
    Rcpp::IntegerMatrix mat_trans = b_make_transition_by_quality_matrix(bb, subs, has_quals && use_quals, err.ncol());
    Rcpp::NumericMatrix mat_quals = b_make_cluster_quality_matrix(bb, subs, has_quals, maxlen);
    Rcpp::DataFrame df_birth_subs = b_make_birth_subs_df(bb, birth_subs, has_quals);

    // Make map from uniques to cluster
    Rcpp::IntegerVector Rmap(nraw);
    for(i=0;i<bb->nclust;i++) {
      for(r=0;r<bb->bi[i]->nraw;r++) {
        raw = bb->bi[i]->raw[r];
        if(raw->correct) { 
          Rmap(raw->index) = i+1; // +1 for R 1-indexing; 
        } else {
          Rmap(raw->index) = NA_INTEGER;
        }
      }
    }
    
    // Organize return List  
    rval = Rcpp::List::create(_["clustering"] = df_clustering, _["birth_subs"] = df_birth_subs, _["subqual"] = mat_trans, _["clusterquals"] = mat_quals, _["map"] = Rmap, _["pval"] = Rpraw);
  } catch(...) {
    cleanup();
    throw;
  }
  cleanup();
  return rval;
}

void run_dada(B *bb, Rcpp::NumericMatrix errMat, 
              int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, 
              bool detect_singletons,
              int max_clust, double min_fold, int min_hamming, int min_abund, 
              bool vectorized_alignment, bool multithread, bool verbose, 
              int SSE, bool gapless, bool greedy, AlignCache *cache, const double *emax) {
  int newi=0, nshuffle = 0;
  bool shuffled = false;

  // Everyone gets aligned within the initial cluster, no KMER screen (and no cache: see AlignCache)
  if(multithread) { b_compare_parallel(bb, 0, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless, greedy, verbose, NULL, NULL); }
  else { b_compare(bb, 0, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless, greedy, verbose, NULL, NULL); }
  b_p_update(bb, greedy, detect_singletons); // Calculates abundance p-value for each raw in its cluster (consensuses)
  
  if(max_clust < 1) { max_clust = bb->nraw; }
  
  while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, min_abund, verbose)) ) {
    if(verbose) Rprintf("\nNew Cluster C%i:", newi);
    if(multithread) { b_compare_parallel(bb, newi, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless, greedy, verbose, cache, emax); }
    else { b_compare(bb, newi, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless, greedy, verbose, cache, emax); }
    // Keep shuffling and updating until no more shuffles
    nshuffle = 0;
    do {
      shuffled = b_shuffle2(bb);
      if(verbose) { Rprintf("S"); }
    } while(shuffled && ++nshuffle < MAX_SHUFFLE);
    if(verbose && nshuffle >= MAX_SHUFFLE) { Rprintf("Warning: Reached maximum (%i) shuffles.\n", MAX_SHUFFLE); }

    b_p_update(bb, greedy, detect_singletons);
    Rcpp::checkUserInterrupt();
  } // while
  
  if(verbose) Rprintf("\nALIGN: %i aligns, %i shrouded (%i raw).\n", bb->nalign, bb->nshroud, bb->nraw);
}

// DadaS: creates an (empty) alignment cache, returned to R as an external pointer. Freed by R's garbage collector.
// [[Rcpp::export]]
SEXP dada_align_cache_new() {
  Rcpp::XPtr<AlignCache> xp(new AlignCache(), true);
  return xp;
}

// DadaS: number of centers and bytes held by an alignment cache (diagnostics).
// [[Rcpp::export]]
Rcpp::NumericVector dada_align_cache_size(SEXP align_cache) {
  Rcpp::XPtr<AlignCache> xp(align_cache);
  double bytes = 0;
  for(auto &kv : xp->by_center) { bytes += kv.second.size(); }
  return Rcpp::NumericVector::create(_["centers"] = (double) xp->by_center.size(), _["bytes"] = bytes);
}
