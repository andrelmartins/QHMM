#include "param_record.hpp"
#include <cassert>

typedef std::vector<double> * vdbl_ptr;

void ParamRecord::init(Params * par) {
  _size = 0;
  _indexes = NULL;
  _records = NULL;

  if (par != NULL) {
    _size = par->length();
    
    for (int i = 0; i < _size; ++i)
      if (par->isFixed(i))
	--_size;
    
    if (_size > 0) {
      _indexes = new int[_size];
      _records = new vdbl_ptr[_size];

      for (int i = 0, j = 0; i < _size; ++i)
	if (!par->isFixed(i))
	  _indexes[j++] = i;

      for (int i = 0; i < _size; ++i)
	_records[i] = new std::vector<double>();
    }

    delete par;
  }
}

ParamRecord::ParamRecord(TransitionFunction * func) : isTransition(true), stateID(func->stateID()), slotID(-1) {
  _funcTrans = func;
  Params * par = func->getParams();
  init(par);
}


ParamRecord::ParamRecord(EmissionFunction * func) : isTransition(false), stateID(func->stateID()), slotID(func->slotID()) {
  _funcEmiss = func;
  Params * par = func->getParams();
  init(par);
}

ParamRecord::~ParamRecord() {
  if (_size > 0) {
    delete[] _indexes;
    
    for (int i = 0; i < _size; ++i)
      delete _records[i];
    delete[] _records;
  }
}

void ParamRecord::collect() {
  if (_size == 0)
    return;

  Params * par;
  if (isTransition)
    par = _funcTrans->getParams();
  else
    par = _funcEmiss->getParams();

  /* add parameters to record */
  for (int i = 0; i < _size; ++i) {
    double val = (*par)[i];
    _records[i]->push_back(val);
  }

  /* */
  delete par;
}

double ParamRecord::value(int row, int position) const {
  assert(row > 0 && row < this->size());
  assert(position >= 0 && position < _size);

  return (*(_records[position]))[row];
}

int ParamRecord::size() const {
  if (_size == 0)
    return 0;
  return _records[0]->size();
}

int ParamRecord::paramSize() const {
  return _size;
}

int ParamRecord::paramIndex(int position) const {
  assert(position >= 0 && position < _size);

  return _indexes[position];
}
