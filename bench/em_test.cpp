#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;

class EM_SS {
public:
  virtual ~EM_SS() {}
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

void test_stream(vector<EM_SS*> & ss_states, int n_states, double ** posts, double * emissions, unsigned int length) {
  
  for (unsigned int i = 0; i < n_states; ++i)
    ss_states[i]->reset();

  for (unsigned int i = 0; i < n_states; ++i)
    ss_states[i]->ss_collect_stream(posts[i], emissions, length);
}

void test_step(vector<EM_SS*> & ss_states, unsigned int n_states, double * posterior, double * emissions, unsigned int length) {

  for (unsigned int i = 0; i < n_states; ++i)
    ss_states[i]->reset();
  
  for (unsigned int i = 0; i < length; ++i, ++emissions) {
    for (unsigned int j = 0; j < n_states; ++j, ++posterior)
      ss_states[j]->ss_collect_step(posterior, emissions);
  }
}

double ** generate_data(unsigned int n_states, unsigned int length, double ** out_table) {
  double ** out;
  out = new double*[n_states];
  *out_table = new double[n_states * length];
  
  for (unsigned int i = 0; i < n_states; ++i) {
    out[i] = new double[length];
    
    for (unsigned int j = 0; j < length; ++j) {
      double u = drand48(); // value in [0, 1)

      out[i][j] = u;
      (*out_table)[n_states * j + i] = u;
    }
    
  }
  
  return out;
}

int main(int argc, char ** argv) {
  unsigned int length;
  unsigned int n_states;
  vector<EM_SS*> ss_states;
  double * emissions;
  double ** posts;
  double * post_tbl;
  double time1, time2;
  
  /* length must be bigger than the combined CPU cache sizes */
//length = 16000000; /* 1600k * 8 ~ 120 MB */
//length = 1000;
  length = atoi(argv[1]);
  n_states = atoi(argv[2]);
  
  /*  EM SS vector */
  for (unsigned int i = 0; i < n_states; ++i) {
    if (i % 2 == 0)
      ss_states.push_back(new EM_SS_A());
    else
      ss_states.push_back(new EM_SS_B());
  }
  
  /* generate data */
  srand48((unsigned long)time(NULL));
  emissions = new double[length];
  
  for (unsigned int i = 0; i < length; ++i)
    emissions[i] = i % 10;
  
  posts = generate_data(n_states, length, &post_tbl);
  
  /* tests */
  clock_t start = clock();
  test_stream(ss_states, n_states, posts, emissions, length);
  clock_t end = clock();

  time1 = (end - start)*1000.0/CLOCKS_PER_SEC;
  
  start = clock();
  test_step(ss_states, n_states, post_tbl, emissions, length);
  end = clock();
  
  time2 = (end - start)*1000.0/CLOCKS_PER_SEC;
  
  cout << "Time(stream)\tTime(step)\tLength" << endl;
  cout << time1 << "\t" << time2 << "\t" << length << endl;
  
  /* clean up */
  delete[] emissions;
  delete[] post_tbl;
  for (unsigned int i = 0; i < n_states; ++i) {
    delete[] posts[i];
    delete ss_states[i];
  }
  delete[] posts;
  
  return EXIT_SUCCESS;
}
