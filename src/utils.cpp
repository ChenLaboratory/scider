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

// Robust double comparison that avoids floating point error based on code by 
// Bruce Dawson.
// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
bool LargerOrAlmostEqual(double A, double B,
                         double maxDiff,
                         double maxRelDiff)
{
  //// A is clearly larger than B.
  if (A > B) return true;
  
  //// A is almost equal to B.
  // Absolute difference. Needed when A,B are close to 0. 
  double diff = fabs(A - B);
  if (diff <= maxDiff)
    return true;
  
  // Relative difference. When A,B are far away from 0. 
  // Should rarely hit here since we z-standardize data before Moran's I.
  A = fabs(A);
  B = fabs(B);
  double largest = (B > A) ? B : A;
  if (diff <= largest * maxRelDiff)
    return true;
  return false;
}