#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;

class EM_SS {
public:
  virtual void reset() = 0;
  
  virtual void ss_collect_stream(double * posterior, double * emissions, unsigned int length) = 0;
  virtual void ss_collect_step(double * posterior_i, double * emissions_i) = 0;
  
  virtual double value() = 0;
};

class EM_SS_A : public EM_SS {
public:
  void reset() {
    _sum_prod = 0;
    _sum = 0;
  }
  
  void ss_collect_stream(double * posterior, double * emissions, unsigned int length) {
    for (unsigned int i = 0; i < length; ++i) {
      double val1 = posterior[i];
      double val2 = emissions[i];
      
      _sum_prod += val1 * val2;
      _sum += val1;
    }
  }
  
  void ss_collect_step(double * posterior_i, double * emissions_i) {
    double val1 = *posterior_i;
    double val2 = *emissions_i;
    
    _sum_prod += val1 * val2;
    _sum += val1;
  }
  
  double value() {
    return _sum_prod / _sum;
  }
  
private:
  double _sum_prod;
  double _sum;
};

class EM_SS_B : public EM_SS {
  void reset() {
    _sum_prod2 = 0;
    _sum = 0;
  }

  void ss_collect_stream(double * posterior, double * emissions, unsigned int length) {
    for (unsigned int i = 0; i < length; ++i) {
      double val1 = posterior[i];
      double val2 = emissions[i];
      
      _sum_prod2 += val1 * val2 * val2;
      _sum += val1;
    }
  }
  
  void ss_collect_step(double * posterior_i, double * emissions_i) {
    double val1 = *posterior_i;
    double val2 = *emissions_i;
    
    _sum_prod2 += val1 * val2 * val2;
    _sum += val1;
  }
  
  double value() {
    return _sum_prod2 / _sum;
  }
  
  
private:
  double _sum_prod2;
  double _sum;
};

void test_stream(EM_SS * ss1, EM_SS * ss2, double * posterior1, double * posterior2, double * emissions, unsigned int length) {
  ss1->reset();
  ss2->reset();
  
  ss1->ss_collect_stream(posterior1, emissions, length);
  ss2->ss_collect_stream(posterior2, emissions, length);
  
  cout << "Stream: " << "[1] = " << ss1->value() << " [2] = " << ss2->value() << endl;
}

void test_step(EM_SS * ss1, EM_SS * ss2, double * posterior12, double * emissions, unsigned int length) {
  ss1->reset();
  ss2->reset();
  
  for (unsigned int i = 0; i < length; ++i) {
    ss1->ss_collect_step(posterior12, emissions);
    ++posterior12;
    ss2->ss_collect_step(posterior12, emissions);
    ++posterior12;
    ++emissions;
  }
  
  cout << "Step: " << "[1] = " << ss1->value() << " [2] = " << ss2->value() << endl;
}

int main(int argc, char ** argv) {
  unsigned int length;
  EM_SS_A ss1;
  EM_SS_B ss2;
  double * emissions;
  double * post1;
  double * post2;
  double * post12;
  double time1, time2;
  
  /* length must be bigger than the combined CPU cache sizes */
//length = 16000000; /* 1600k * 8 ~ 120 MB */
//length = 1000;
  length = atoi(argv[1]);
  
  /* generate data */
  emissions = new double[length];
  post1 = new double[length];
  post2 = new double[length];
  post12 = new double[2*length];
  
  srand48((unsigned long)time(NULL));

  for (unsigned int i = 0; i < length; ++i) {
    double u = drand48(); // value in [0, 1)

    post1[i] = u;
    post2[i] = 1.0 - u;
    post12[2*i] = u;
    post12[2*i + 1] = 1.0 - u;
    emissions[i] = i % 10;
  }
  
  /* tests */
  clock_t start = clock();
  test_stream(&ss1, &ss2, post1, post2, emissions, length);
  clock_t end = clock();

  time1 = (end - start)*1000.0/CLOCKS_PER_SEC;
  
  start = clock();
  test_step(&ss1, &ss2, post12, emissions, length);
  end = clock();
  
  time2 = (end - start)*1000.0/CLOCKS_PER_SEC;
  
  cout << "Time(stream)\tTime(step)\tLength" << endl;
  cout << time1 << "\t" << time2 << "\t" << length << endl;
  
  /* clean up */
  delete[] emissions;
  delete[] post1;
  delete[] post2;
  delete[] post12;
  
  return EXIT_SUCCESS;
}
