#define R_NO_REMAP
#include "utils.h"

void pseudoP(int start,int end,uint64_t seed_start,
             R_xlen_t n,
             int *p_n_nbrs,
             double **weight,
             int perms,
             int max_rand,
             double *p_data2,
             double *p_data1,
             double *lisa_vec,
             double *sig_local_vec,
             double *cluster_vec,
             double cutoff) {
  // possible random neighbors must has at least 1 neighbors
  int n_sample = 0; 
  int *potential_nbrs = new int[n]; //size n but will only use the first n_sample values.
  int seed_offset = 0;
  for (int i=0;i<n;i++) {
    if (p_n_nbrs[i]>0) {
      potential_nbrs[n_sample]=i;
      n_sample++;
    }
    seed_offset += (i<start)*p_n_nbrs[i];
  }
  seed_start += seed_offset*perms;
  // index for sampling. Each point cannot be its own neighbors hence, n_sample--
  n_sample--;
  int *x = new int[n_sample];
  for (int cnt = start; cnt <= end;cnt++) {
    if (cluster_vec[cnt] == 6 || cluster_vec[cnt] == 5) {// Skip CLUSTER_NEIGHBORLESS||CLUSTER_UNDEF
      continue;
    }
    
    for (int i=0;i<n_sample;i++) x[i]=i;
    uint64_t countLarger = 0;
    
    for(int p=0; p<perms; p++) {
      double permuted_lag = 0;
      // Calculate moran for random nbrs
      sample_without_replacement(x,n_sample,p_n_nbrs[cnt],seed_start);
      seed_start += p_n_nbrs[cnt];
      for (int i=n_sample-p_n_nbrs[cnt]; i<=n_sample-1; i++) {
        if (potential_nbrs[x[i]]==cnt) x[i] = n_sample;
        permuted_lag += p_data2[potential_nbrs[x[i]]]*weight[cnt][i-(n_sample-p_n_nbrs[cnt])];
      }
      countLarger += (permuted_lag*p_data1[cnt] > lisa_vec[cnt]);
    }
    if (perms-countLarger <= countLarger) {
      countLarger = perms-countLarger;
    }
    
    sig_local_vec[cnt] = (countLarger+1.0)/(perms+1);

    if (sig_local_vec[cnt] > cutoff) {
      cluster_vec[cnt] = 0; // CLUSTER_NOT_SIG
    }
  }
  delete[] x;
  delete[] potential_nbrs;
}

extern "C" {
  SEXP C_localMoran(SEXP nbrs,
                    SEXP n_nbrs,
                    SEXP weight_,
                    SEXP data1,
                    SEXP data2,
                    SEXP significance_cutoff,
                    SEXP permutations,
                    SEXP seed,
                    SEXP cpu_threads,
                    SEXP hhonly_) {
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
    double cutoff = Rf_asReal(significance_cutoff);
    int perms = Rf_asInteger(permutations);
    int n_cpu = Rf_asInteger(cpu_threads);
    bool hhonly = Rf_asLogical(hhonly_);
    
    // Standardize data
    p_data1 = standardizeData(p_data1,n);
    p_data2 = standardizeData(p_data2,n);
    
    // Calculating Moran
    double* lag_vec = new double[n]();
    double* lisa_vec = new double[n];
    
    // int* nbrs_i;
    for (R_xlen_t i=0; i < n; i++) {
      if (p_n_nbrs[i] == 0) {
        continue;
      }
      for (int j = 0; j < p_n_nbrs[i];j++) {
        lag_vec[i] += p_data2[p_nbrs[i][j]-1]*weight[i][j];
      }
      lisa_vec[i] = lag_vec[i]*p_data1[i];
    }
    
    // Assigning clusters
    double* cluster = new double[n];
    for (R_xlen_t i=0; i < n; i++) {
      if(p_n_nbrs[i] > 0) {
        if (p_data1[i] > 0  && lag_vec[i] < 0) {cluster[i] = 4;} // CLUSTER_HIGHLOW
        else if (p_data1[i] < 0  && lag_vec[i] > 0) {cluster[i] = 3;} // CLUSTER_LOWHIGH
        else if (p_data1[i] < 0  && lag_vec[i] < 0) {cluster[i] = 2;} // CLUSTER_LOWLOW
        else {cluster[i] = 1;} // CLUSTER_HIGHHIGH
        
        if (hhonly && cluster[i]!=1) {
          cluster[i] = 5; // CLUSTER_UNDEF
        }
      } else {
        cluster[i] = 6; // CLUSTER_NEIGHBORLESS
      }
    }
    
    // Multithreading pseudo p-value
    double* sig_local_vec = new double[n](); // store p-value output
    int max_rand = n-1;
    parallel(n_cpu,n,[&](int start, int end) -> void {
      pseudoP(start,end,
              p_seed,
              n,
              p_n_nbrs,
              weight,
              perms,
              max_rand,
              p_data2,
              p_data1,
              lisa_vec,
              sig_local_vec,
              cluster,
              cutoff);
    });
    
    // Wrangling C to R
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    
    SEXP names = PROTECT(Rf_allocVector(STRSXP,4));
    SET_STRING_ELT(names,0,Rf_mkChar("lisa"));
    SET_STRING_ELT(names,1,Rf_mkChar("cluster"));
    SET_STRING_ELT(names,2,Rf_mkChar("p.value"));
    SET_STRING_ELT(names,3,Rf_mkChar("spatiallag"));
    Rf_setAttrib(out, Rf_install("names"), names);
    
    SEXP lisa_vec_final = SET_VECTOR_ELT(out,0,Rf_allocVector(REALSXP,n));
    SEXP cluster_vec_final = SET_VECTOR_ELT(out,1,Rf_allocVector(INTSXP,n));
    SEXP sig_local_vec_final = SET_VECTOR_ELT(out,2,Rf_allocVector(REALSXP,n));
    SEXP lag_vec_final = SET_VECTOR_ELT(out,3,Rf_allocVector(REALSXP,n));
    
    memcpy(REAL(lag_vec_final),lag_vec,n*sizeof(double));
    
    double* lisa_vec_final_p = REAL(lisa_vec_final);
    int* cluster_vec_final_p = INTEGER(cluster_vec_final);
    double* sig_local_vec_final_p = REAL(sig_local_vec_final);
    
    for (R_xlen_t i=0; i<n; i++){
      lisa_vec_final_p[i]=lisa_vec[i];
      cluster_vec_final_p[i]=cluster[i];
      sig_local_vec_final_p[i]=sig_local_vec[i];
    }
    
    delete [] p_nbrs;
    delete [] weight;
    delete [] p_data1;
    delete [] p_data2;
    delete [] lag_vec;
    delete [] lisa_vec;
    delete [] cluster;
    delete [] sig_local_vec;
    
    UNPROTECT(7);
    return(out);
  }
}