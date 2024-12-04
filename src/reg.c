#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP C_localMoran(SEXP nbrs,
                  SEXP n_nbrs,
                  SEXP data1,
                  SEXP data2,
                  SEXP significance_cutoff,
                  SEXP permutations,
                  SEXP seed);
static const R_CallMethodDef CallEntries[] = {
  {"C_localMoran", (DL_FUNC) &C_localMoran,7},
  {NULL, NULL, 0}
};


void R_init_hexDensity(DllInfo *info) {
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
