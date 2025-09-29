#include "utils.h"

// standardize data
double* standardizeData(double* data, R_xlen_t n) {
  double* res = new double[n];
  double mean = 0;
  for (R_xlen_t i=0;i<n;i++) {
    mean += data[i];
  }
  mean /= n;
  for (R_xlen_t i=0;i<n;i++) {
    res[i] = data[i] - mean;
  }
  double sd = 0.0;
  for (R_xlen_t i=0;i<n;i++) {
    sd += res[i]*res[i];
  }
  sd = std::sqrt(sd/(double)(n-1.0));
  for (R_xlen_t i=0;i<n;i++) {
    res[i] /= sd;
  }
  return(res);
}

// from GeoDaCenter/libgeoda
// TODO: Thomas Wang is apparently not that great??
double ThomasWangHashDouble(uint64_t key) {
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return 5.42101086242752217E-20 * key;
}

// This is a weird implementation that is highly specific to the task.
// Guarantees sampling k values within k picks (Important for stable multi-threading since seed will be predictable).
// Samples will be stored in the range x[n-k]...x[n-1] inclusive.
// x will be (partially) permuted after each run.
// template <typename T> void sample_without_replacement(T* x, int n,int k, uint64_t seed)
void sample_without_replacement(int* x,
                                int n,
                                int k,
                                uint64_t seed)
  {
  int tempSwap;
  int j;
  double rng_val;
  for (int i=0;i<k;i++) {
    rng_val = ThomasWangHashDouble(seed++)*(--n);
    j = (int)(rng_val<0.0?ceil(rng_val - 0.5):floor(rng_val + 0.5));
    
    // Swap
    tempSwap = x[j];
    x[j] = x[n];
    x[n] = tempSwap;
  }
}

// n choose k (w/out replacement) from range 0:n
int* sample_to_n(int n,
                 int k,
                 uint64_t seed)
{
  int* res = new int[k];
  int* chosen = new int[n](); // if has chosen then go choose this one instead.  
  int j;
  double rng_val;
  for (int i=0;i<k;i++) {
    rng_val = ThomasWangHashDouble(seed++)*(--n);
    j = (int)(rng_val<0.0?ceil(rng_val - 0.5):floor(rng_val + 0.5));
    
    res[i] = chosen[j]?chosen[j]:j;
    chosen[j] = n;
  }
  
  delete [] chosen;
  return(res);
}

// Break n_task into appropriate chunks according to n_cpu then pass the 'start'
// and 'end' of each chunk (inclusive) to f.
void parallel(int n_cpu,int n_task, std::function<void(int, int)> f) {
  // Getting number of threads
  int quotient = n_task/n_cpu;
  int remainder = n_task%n_cpu;
  int tot_threads = quotient?n_cpu:remainder;
  // Parallel
  std::thread *threads = new std::thread[tot_threads];
  
  for (int i=0; i<tot_threads; i++) {
    int a,b;
    if (i < remainder) {
      a = i*(quotient+1);
      b = a+quotient;
    } else {
      a = i*quotient+remainder;//remainder*(quotient+1) + (i-remainder)*quotient;
      b = a+quotient-1;
    }
    threads[i] = std::thread(f,a,b);
  }
  for (int i=0;i<tot_threads;i++) { 
    threads[i].join();
  }
  delete[] threads;
}