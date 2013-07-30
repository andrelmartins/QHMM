#include <iostream>
#include <cstdlib>
#include <ctime>
#include "../src/logsum.hpp"

using namespace std;

void add_numbers(LogSum * lg, int N, double fraction) {
  for (int i = 0; i < N; ++i) {
    double u = drand48(); // value in [0, 1)
    // on average "fraction" % of the numbers are lower than threshold    
    double log_value = LogSum::SUM_LOG_THRESHOLD / (1.0 - fraction) * u;
    
    lg->store(log_value);
  }
}

void run_test(LogSum * lg, int repeats, int N, double fraction, const char * name) {
  double result = 0; 
  
  clock_t start = clock();
  for (int k = 0; k < repeats; ++k) {
    lg->clear();
    
    add_numbers(lg, N, fraction);
    
    result += lg->compute();
  }
  clock_t end = clock();
  
  cout << name << ": " << result << " in " << (end - start)*1000.0/CLOCKS_PER_SEC << " ms" << endl;
}

int main(int argc, char ** argv) {
  int size;
  int repeats;
  double fraction;
  
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " <size> <repeats> <drop fraction>\n" << endl;
    return EXIT_FAILURE;
  }
  
  size = atoi(argv[1]);
  repeats = atoi(argv[2]);
  fraction = atof(argv[3]);
  
  /* output time for each implementation type */
  srand48((unsigned long)time(NULL));
  
  /* optimized and non-optimized version */
  
  if (size == 2) { /* type 1 */
    {
      LogSum * lg2 = LogSum::createType(1, 2, true);
      run_test(lg2, repeats, 2, fraction, "LogSum2(opt)");
      delete lg2;
    }
    {
      LogSum * lg2 = LogSum::createType(1, 2, false);
      run_test(lg2, repeats, 2, fraction, "LogSum2");
      delete lg2;
    }
  }
  
  /* optimized */
  {
    LogSum * lg = LogSum::createType(2, size, true);
    run_test(lg, repeats, size, fraction, "LogSum(opt)");
    delete lg;
  }
  
  /* non-optimized */
  {
    LogSum * lg = LogSum::createType(2, size, false);
    run_test(lg, repeats, size, fraction, "LogSum");
    delete lg;
  }
  
  return EXIT_SUCCESS;
}
