#include "utils.h"
#include "rand.h"

// Bootstrap Moran's I by shuffling all points around
// PCG-PRNG
void pseudoP_global(int start, int end, uint64_t seed_start,
                    R_xlen_t n,
                    double W,
                    int *p_n_nbrs,
                    double **weight,
                    int **p_nbrs,
                    double *p_data2,
                    double *p_data1,
                    double* permuted_moran){
  // Initiate rng & figure out correct starting position for the thread
  pcg32_random_t rng;
  pcg32_srandom_r(&rng, seed_start);
  rng.state = pcg_advance_lcg_64(rng.state,start*n);
  
  // Start bootstrap
  int *x = new int[n];
  for (int p=start;p<=end;p++) {
    // Permuting all points (n chooses n).
    for (R_xlen_t i=0;i<n;i++) x[i]=i;
    sample_without_replacement(x,n,n,&rng);

    // Calculating Moran's I
    for (R_xlen_t i=0; i < n; i++) {
      if (p_n_nbrs[i] == 0) {
        continue;
      }
      double lag = 0;
      for (int j = 0; j < p_n_nbrs[i];j++) {
        lag += p_data2[x[p_nbrs[i][j]-1]]*weight[i][j];
      }
      permuted_moran[p] += lag*p_data1[x[i]];
    }
    permuted_moran[p] /= W;
  }
  delete [] x;
}

extern "C" {
  SEXP C_globalMoran(SEXP nbrs,
                     SEXP n_nbrs,
                     SEXP weight_,
                     SEXP data1,
                     SEXP data2,
                     SEXP permutations,
                     SEXP seed,
                     SEXP cpu_threads) {
    // Wrangling R to C
    nbrs = PROTECT(Rf_coerceVector(nbrs,VECSXP));
    n_nbrs = PROTECT(Rf_coerceVector(n_nbrs,INTSXP));
    weight_ = PROTECT(Rf_coerceVector(weight_,VECSXP));
    data1 = PROTECT(Rf_coerceVector(data1,REALSXP));
    data2 = PROTECT(Rf_coerceVector(data2,REALSXP));
    
    R_xlen_t n = Rf_xlength(n_nbrs);
    
    int** p_nbrs = new int*[n];
    double** weight = new double*[n];
    for (R_xlen_t i=0; i<n; i++){
      p_nbrs[i] = INTEGER(Rf_coerceVector((VECTOR_ELT(nbrs,i)),INTSXP));
      weight[i] = REAL(Rf_coerceVector((VECTOR_ELT(weight_,i)),REALSXP));
    } 
    double* p_data1 = (double*) REAL(data1);
    double* p_data2 = (double*) REAL(data2);
    int* p_n_nbrs = INTEGER(n_nbrs);
    
    uint64_t  p_seed = (uint64_t )Rf_asInteger(seed);
    
    int perms = Rf_asInteger(permutations);
    int n_cpu = Rf_asInteger(cpu_threads);
    
    // Shuffle so sequential seeds are different. Not strictly needed.
    // Numbers taken from r-source/src/main/RNG.c
    for(int i = 0; i < 50; i++) p_seed = (69069 * p_seed + 1);
    
    // Standardize data
    p_data1 = standardizeData(p_data1,n);
    p_data2 = standardizeData(p_data2,n);
    
    // Calculating Moran
    double lisa = 0;
    double W = 0;
    for (R_xlen_t i=0; i < n; i++) {
      double lag = 0;
      for (int j = 0; j < p_n_nbrs[i];j++) {
        lag += p_data2[p_nbrs[i][j]-1]*weight[i][j];
        W += weight[i][j];
      }
      lisa += lag*p_data1[i];
    }
    lisa /= W;
    
    // Bootstrap Moran's I
    double* permuted_moran = new double[perms](); // store permuted moran output
    parallel(n_cpu,perms,[&](int start, int end) -> void {
      pseudoP_global(start,end,
                      p_seed,
                      n,
                      W,
                      p_n_nbrs,
                      weight,
                      p_nbrs,
                      p_data2,
                      p_data1,
                      permuted_moran);
    });
    
    // Folded p-value
    int countLarger = 0;
    for (int p=0;p<perms;p++){
      countLarger += permuted_moran[p]>=lisa;
    }
    if (perms-countLarger <= countLarger) {
      countLarger = perms-countLarger;
    }
    double pvalue = (countLarger+1.0)/(perms+1);
    
    // Wrangling C to R
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    
    SEXP names = PROTECT(Rf_allocVector(STRSXP,3));
    SET_STRING_ELT(names,0,Rf_mkChar("lisa"));
    SET_STRING_ELT(names,1,Rf_mkChar("p.value"));
    SET_STRING_ELT(names,2,Rf_mkChar("permuted_lisa"));
    Rf_setAttrib(out, Rf_install("names"), names);
    
    SET_VECTOR_ELT(out,0,Rf_ScalarReal(lisa));
    SET_VECTOR_ELT(out,1,Rf_ScalarReal(pvalue));
    SEXP permuted_lisa_final = SET_VECTOR_ELT(out,2,Rf_allocVector(REALSXP,perms));
    double* permuted_lisa_final_p = REAL(permuted_lisa_final);
    for (R_xlen_t p=0; p<perms; p++){
      permuted_lisa_final_p[p]=permuted_moran[p];
    }
    
    delete [] p_nbrs;
    delete [] weight;
    delete [] p_data1;
    delete [] p_data2;
    delete [] permuted_moran;
    
    UNPROTECT(7);
    return(out);
  }
}