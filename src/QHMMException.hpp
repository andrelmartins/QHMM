#ifndef QHMMEXCEPTION_HPP
#define QHMMEXCEPTION_HPP

#include <exception>
#include <vector>
#include <string>

using namespace std;

class QHMMException : public exception {

public:
  QHMMException(const QHMMException & other) {
    this->is_transition = other.is_transition;
    this->state = other.state;
    this->slot = other.slot;
    this->sequence_index = other.index;
    this->evalue = other.evalue;
    this->sequence_id = other.sequence_id;
    this->stack = other.stack;
    this->_msg = other._msg;
  }

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

#endif
