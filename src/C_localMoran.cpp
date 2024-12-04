#define R_NO_REMAP
#include <Rinternals.h>
#include <math.h> //for sqrt
#include <cstring> //for memset
// from GeoDa
// TODO: Check copyright https://github.com/GeoDaCenter/libgeoda/blob/a2a40be473dd27f7ed60b5c0e9f69e0302ab8d4e/GenUtils.cpp#L513 
// TODO: Thomas Wang is apparently not that great??
double ThomasWangHashDouble(unsigned __int64 key) {
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return 5.42101086242752217E-20 * key;
}

extern "C" {
  SEXP C_localMoran(SEXP nbrs,
                    SEXP n_nbrs,
                    SEXP data1,
                    SEXP data2,
                    SEXP significance_cutoff,
                    SEXP permutations,
                    SEXP seed) {
    nbrs = PROTECT(Rf_coerceVector(nbrs,VECSXP));
    n_nbrs = PROTECT(Rf_coerceVector(n_nbrs,INTSXP));
    data1 = PROTECT(Rf_coerceVector(data1,REALSXP));
    data2 = PROTECT(Rf_coerceVector(data2,REALSXP));
    seed = PROTECT(Rf_coerceVector(seed,INTSXP));
    
    double *p_data1 = (double*) REAL(data1);
    double *p_data2 = (double*) REAL(data2);
    int *p_n_nbrs = INTEGER(n_nbrs);

    unsigned __int64  p_seed = (unsigned __int64 )Rf_asInteger(seed);
    double cutoff = Rf_asReal(significance_cutoff);
    int perms = Rf_asInteger(permutations);
    
    R_xlen_t n = Rf_xlength(n_nbrs);
    
    // Standardize data
    double mean1 = 0.0;
    double mean2 = 0.0;
    for (R_xlen_t i=0; i<n; i++) {
      mean1 += p_data1[i];
      mean2 += p_data2[i];
    }
    mean1 = mean1/n;
    mean2 = mean2/n;
    for (R_xlen_t i=0; i<n; i++) {
      p_data1[i] -= mean1;
      p_data2[i] -= mean2;
    }
    double sd1 = 0.0;
    double sd2 = 0.0;
    for (R_xlen_t i=0; i<n; i++) {
      sd1 += p_data1[i] * p_data1[i];
      sd2 += p_data2[i] * p_data2[i];
    }
    sd1 = std::sqrt(sd1/(double)(n-1.0));
    sd2 = std::sqrt(sd2/(double)(n-1.0));
    for (R_xlen_t i=0; i<n; i++) {
      p_data1[i] /= sd1;
      p_data2[i] /= sd2;
    }

    // Calculating moran
    double* lag_vec = new double[n]();
    double* lisa_vec = new double[n];
    double* cluster_vec = new double[n];
    
    int* nbrs_i;
    for (R_xlen_t i=0; i < n; i++) {
      if (p_n_nbrs[i] == 0) {
        continue;
      }
      nbrs_i = INTEGER(Rf_coerceVector((VECTOR_ELT(nbrs,i)),INTSXP));
      for (int j = 0; j < p_n_nbrs[i];j++) {
        lag_vec[i] += p_data2[nbrs_i[j]-1];
      }
      lag_vec[i] /= p_n_nbrs[i];
      lisa_vec[i] = lag_vec[i]*p_data1[i];
    }

    // Assigning clusters
    for (R_xlen_t i=0; i < n; i++) {
      if(p_n_nbrs[i] > 0) {
        if (p_data1[i] > 0  && lag_vec[i] < 0) {cluster_vec[i] = 4;} // CLUSTER_HIGHLOW
        else if (p_data1[i] < 0  && lag_vec[i] > 0) {cluster_vec[i] = 3;} // CLUSTER_LOWHIGH
        else if (p_data1[i] < 0  && lag_vec[i] < 0) {cluster_vec[i] = 2;} // CLUSTER_LOWLOW
        else {cluster_vec[i] = 1;} // CLUSTER_HIGHHIGH
      } else {
        cluster_vec[i] = 6; // CLUSTER_NEIGHBORLESS
      }
    }
    
    // Pseudo p-value
    double* sig_local_vec = new double[n]();
    int max_rand = n-1;
    
    // Set for random nbrs
    int current = 0;
    int* buffer = new int [n];
    char* flags = new char [n];
    std::memset(flags, '\x0', n);
    
    double* permutedSA = new double[n]();
    unsigned __int64 countLarger;
    for (R_xlen_t cnt = 0; cnt < n;cnt++) {
      if (p_n_nbrs[cnt] == 0) {
        continue;
      }

      for(int p=0; p<perms; p++) {
        // Get random nbrs
        int rand=0, newRandom;
        double rng_val;
        while (rand < p_n_nbrs[cnt]) {
          rng_val = ThomasWangHashDouble(p_seed++) * max_rand;
          newRandom = (int)(rng_val<0.0?ceil(rng_val - 0.5):floor(rng_val + 0.5));
          
          if (newRandom!=cnt && (flags[newRandom]==0) && p_n_nbrs[newRandom]>0) {
            // set push
            buffer[current++] = newRandom;
            flags[newRandom] = 'i'; //TODO: any reason for "i" specifically???
            rand++;
          }
        }
        
        double permuted_lag = 0;
        // Calculate moran for random nbrs
        for (int n=0; n<p_n_nbrs[cnt]; n++) {
          // set pop
          permuted_lag += p_data2[buffer[--current]];
          flags[buffer[current]] = '\x0';
        }
        
        if (p_n_nbrs[cnt]>0) { 
          permuted_lag /= p_n_nbrs[cnt];
        }
        permutedSA[p] = permuted_lag*p_data1[cnt];
      }
      
      countLarger = 0;
      for (int p=0; p<perms; p++) {
        if (permutedSA[p] >= lisa_vec[cnt]) {
          countLarger += 1;
        }
      }
      if (perms-countLarger <= countLarger) {
        countLarger = perms-countLarger;
      }
      sig_local_vec[cnt] = (countLarger+1.0)/(perms+1);
    }
    // Compare p-value for significance
    for (R_xlen_t i=0; i < n; i++) {
      if (sig_local_vec[i] > cutoff ) {
        cluster_vec[i] = 0; // CLUSTER_NOT_SIG
      }
    }

    
    // Formatting output to return to R
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    
    SEXP names = PROTECT(Rf_allocVector(STRSXP,3));
    SET_STRING_ELT(names,0,Rf_mkChar("lisa_vec"));
    SET_STRING_ELT(names,1,Rf_mkChar("cluster_vec"));
    SET_STRING_ELT(names,2,Rf_mkChar("local_pseudo_p"));
    Rf_setAttrib(out, Rf_install("names"), names);
    
    SEXP lisa_vec_final = SET_VECTOR_ELT(out,0,Rf_allocVector(REALSXP,n));
    SEXP cluster_vec_final = SET_VECTOR_ELT(out,1,Rf_allocVector(INTSXP,n));
    SEXP sig_local_vec_final = SET_VECTOR_ELT(out,2,Rf_allocVector(REALSXP,n));
    double* lisa_vec_final_p = REAL(lisa_vec_final);
    int* cluster_vec_final_p = INTEGER(cluster_vec_final);
    double* sig_local_vec_final_p = REAL(sig_local_vec_final);
    
    for (R_xlen_t i=0; i<n; i++){
      lisa_vec_final_p[i]=lisa_vec[i];
      cluster_vec_final_p[i]=cluster_vec[i];
      sig_local_vec_final_p[i]=sig_local_vec[i];
    }
    
    UNPROTECT(7);
    return(out);
  }
}