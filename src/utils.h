#ifndef UTILS_H
#define UTILS_H

#define R_NO_REMAP

#include <Rinternals.h> //for R_xlen_t
#include <cmath> //for sqrt
#include <cstdint> //for uint64_t
#include <thread>
#include <functional>
#include <cstring> // for memcpy
#include <float.h> // for DBL_EPSILON

double* standardizeData(double* data, R_xlen_t n);
void parallel(int n_cpu,int n_task, std::function<void(int, int)> f);
bool LargerOrAlmostEqual(double A, double B, double maxDiff = std::sqrt(DBL_EPSILON), 
                         double maxRelDiff = std::sqrt(DBL_EPSILON));

#endif