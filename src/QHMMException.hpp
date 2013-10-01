#ifndef QHMMEXCEPTION_HPP
#define QHMMEXCEPTION_HPP

#include <exception>
#include <vector>
#include <string>
#include <omp.h>

using namespace std;

class QHMMException : public exception {

public:
  QHMMException(string msg, string frame, bool is_transition, int state, int slot, int index, double value) : _msg(msg) {
    this->is_transition = is_transition;
    this->state = state;
    this->slot = slot;
    this->sequence_index = index;
    this->evalue = value;

    this->sequence_id = -1;
    this->stack.push_back(frame);
  }
  virtual ~QHMMException() throw() {}

  // overriden what() method from exception class
  const char* what() const throw() { return _msg.c_str(); }

  // exception properties
  bool is_transition;
  int state;
  int slot;
  double evalue;
  int sequence_index;
  int sequence_id;

  vector<string> stack;

private:
  string _msg;
};

/*
class ThreadException {
  exception_ptr Ptr;
  mutex         Lock;
public:
  ThreadException(): Ptr(nullptr) {}
  ~ThreadException(){ this->Rethrow(); }  
  void Rethrow(){
    if(this->Ptr) rethrow_exception(this->Ptr);
  }
  void CaptureException() { 
    unique_lock<mutex> guard(this->Lock);
    this->Ptr = std::current_exception(); 
  }
  template <typename Function, typename... Parameters>
  void Run(Function f, Parameters... params)
  {
    try 
    {
      f(params...);
    }
    catch (...)
    {
        CaptureException();
    }
  }

};
*/

#endif
