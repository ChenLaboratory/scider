#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP C_localMoran(SEXP nbrs,
                  SEXP n_nbrs,
                  SEXP weight_,
                  SEXP data1,
                  SEXP data2,
                  SEXP significance_cutoff,
                  SEXP permutations,
                  SEXP seed,
                  SEXP cpu_threads,
                  SEXP hhonly_);
SEXP C_globalMoran(SEXP nbrs,
                   SEXP n_nbrs,
                   SEXP weight_,
                   SEXP data1,
                   SEXP data2,
                   SEXP permutations,
                   SEXP seed,
                   SEXP cpu_threads);
SEXP C_findSNN(SEXP knn,
               SEXP k_,
               SEXP n_cells_,
               SEXP type_,
               SEXP cpu_threads_);
SEXP C_corDensityPerm(SEXP prob,
                      SEXP obs_,
                      SEXP acc_,
                      SEXP b_start_,
                      SEXP nb_,
                      SEXP seed_,
                      SEXP cpu_threads);

static const R_CallMethodDef CallEntries[] = {
    {"C_localMoran", (DL_FUNC) &C_localMoran,10},
    {"C_globalMoran", (DL_FUNC) &C_globalMoran,8},
    {"C_findSNN", (DL_FUNC) &C_findSNN,5},
    {"C_corDensityPerm", (DL_FUNC) &C_corDensityPerm,7},
    {NULL, NULL, 0}
};

void R_init(DllInfo *info) {
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}