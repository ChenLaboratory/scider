#include <vector>
#include <string>
#include "utils.h"
// Contruct Shared Nearest-Neighbor graph with possible weights being rank, 
// jaccard, and number.
// Based on Aaron Lun's code. See https://github.com/libscran/scran_graph_cluster
extern "C" {
  SEXP C_findSNN(SEXP knn,
                 SEXP k_,
                 SEXP n_cells_,
                 SEXP type_,
                 SEXP cpu_threads_) {
    knn = PROTECT(Rf_coerceVector(knn,INTSXP));
    int *knn_p = INTEGER(knn);
    
    int k = Rf_asInteger(k_);
    int n_cells = Rf_asInteger(n_cells_);
    std::string type = CHAR(Rf_asChar(type_));
    int cpu_threads = Rf_asInteger(cpu_threads_);
    
    // Convert to adjacency list
    std::vector<std::vector<int>> simple_host;
    std::vector<std::vector<std::pair<int,int>>> ranked_host;
    if(type=="rank"){
      ranked_host.resize(n_cells);
      for (int i = 0;i < n_cells;i++) {
        ranked_host[i].emplace_back(i,0); //1st nbrs is itself
        int rank = 1;
        for (int j=0;j<k;j++) {
          ranked_host[knn_p[k*i+j]].emplace_back(i,rank);
          rank++;
        }
      }
    } else {
      simple_host.resize(n_cells);
      for (int i = 0;i < n_cells;i++) {
        simple_host[i].emplace_back(i);
        for (int j = 0; j < k; j++) {
          // TODO: make a get_neighbor helper function.
          simple_host[knn_p[k*i+j]].emplace_back(i); //1st nbrs is itself
        }
      } 
    }
    
    // Construct SNN graph in parallel
    std::vector<std::vector<int>> edge_stores(n_cells);
    std::vector<std::vector<double>> weight_stores(n_cells);
    parallel(cpu_threads,n_cells,[&](int start, int end) -> void {
      std::vector<int> current_score(n_cells);
      std::vector<int> current_added;
      current_added.reserve(n_cells); // Skip zero-ing vector + initial size = 0
      
      for (int i = start; i <= end;i++) {
        for (int j = 0; j <= k; j++) {
          int cur_nbrs = (j==0 ? i : knn_p[k*i+j-1]); //1st nbrs is itself
          
          if (type=="rank") {
            for (const auto& host : ranked_host[cur_nbrs]) {
              int othernode = host.first;
              if (othernode != i) {// avoid storing self as nbrs
                int& existing_other = current_score[othernode];
                // Record lowest-combined rank per nbrs
                int currank = host.second + j;
                if (existing_other==0) {
                  existing_other=currank;
                  current_added.emplace_back(othernode);
                } else if (existing_other>currank) {
                  existing_other=currank;
                }
              }
            }
          } else {
            for (const int& othernode : simple_host[cur_nbrs]) {
              if (othernode < i) {// avoid duplicates from SNN symmetry
                int& existing_other = current_score[othernode];
                // Record number of shared neighbors with i
                if (existing_other == 0) {
                  current_added.emplace_back(othernode);
                }
                existing_other++;
              }
            }
          }
        }
        
        // Store the edges & weights.
        edge_stores[i].reserve(current_added.size());
        weight_stores[i].reserve(current_added.size());
        for (int othernode:current_added) {
          // Weight
          int& otherscore = current_score[othernode];
          double finalscore = (double) otherscore;
          if (type=="rank") {
            finalscore = k - 0.5*finalscore;
          }else if (type=="jaccard") {
            finalscore = finalscore/(2*(k+1)-finalscore); //k+1 because need to include the node itself. TODO: check this if change nn algorithm
          }
          weight_stores[i].emplace_back(std::max(finalscore, 1e-6)); // Ensure positive edge is always recorded
          otherscore = 0;
          // Edge
          edge_stores[i].emplace_back(othernode);
        }
        current_added.clear();
      }
    });
    
    // Formatting output to return to R
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP names = PROTECT(Rf_allocVector(STRSXP,2));
    SET_STRING_ELT(names,0,Rf_mkChar("index"));
    SET_STRING_ELT(names,1,Rf_mkChar("weight"));
    Rf_setAttrib(out, Rf_install("names"), names);
    
    SEXP edges_final = SET_VECTOR_ELT(out,0,Rf_allocVector(VECSXP,n_cells));
    SEXP weights_final = SET_VECTOR_ELT(out,1,Rf_allocVector(VECSXP,n_cells));
    
    for (int i=0;i<n_cells;i++) {
      size_t n = edge_stores[i].size();
      SEXP edges_i = SET_VECTOR_ELT(edges_final,i,Rf_allocVector(INTSXP,n));
      SEXP weights_i = SET_VECTOR_ELT(weights_final,i,Rf_allocVector(REALSXP,n));

      memcpy(INTEGER(edges_i),edge_stores[i].data(),n*sizeof(int));
      memcpy(REAL(weights_i),weight_stores[i].data(),n*sizeof(double));
    }
    
    UNPROTECT(3);
    return(out);
  }
}