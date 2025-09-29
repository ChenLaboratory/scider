#define R_NO_REMAP

#include <Rinternals.h> //for R_xlen_t
#include <math.h> //for sqrt
#include <cstdint> //for uint64_t
#include <thread>
#include <functional>
#include <cstring> // for memcpy
double* standardizeData(double* data, R_xlen_t n);
double ThomasWangHashDouble(uint64_t key);
void sample_without_replacement(int* x,int n,int k,uint64_t seed);
int* sample_to_n(int n,int k,uint64_t seed);

void parallel(int n_cpu,int n_task, std::function<void(int, int)> f);