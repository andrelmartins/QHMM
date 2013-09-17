#ifndef QHMMEXCEPTION_HPP
#define QHMMEXCEPTION_HPP

#include <exception>
#include <vector>
#include <string>

using namespace std;

class QHMMException : public exception {

public:
  QHMMException(string msg, bool is_transition, int state, int slot, int index, double value) : _msg(msg) {
    this->is_transition = is_transition;
    this->state = state;
    this->slot = slot;
    this->sequence_index = index;
    this->evalue = value;

    this->sequence_id = -1;
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

#endif
