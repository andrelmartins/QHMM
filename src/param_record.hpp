#ifndef PARAM_RECORD_HPP
#define PARAM_RECORD_HPP

#include <vector>
#include "base_classes.hpp"

class ParamRecord {
public:
  ParamRecord(TransitionFunction * func);
  ParamRecord(EmissionFunction * func);
  ~ParamRecord();

  void collect();
  
  const bool isTransition;
  const int stateID;
  const int slotID;
  
  int size() const;
  int paramSize() const;
  int paramIndex(int position) const;

  double value(int row, int position) const;

private:
  int _size;
  int * _indexes;

  TransitionFunction * _funcTrans;
  EmissionFunction * _funcEmiss;

  std::vector<double> ** _records;

  void init(Params * par);
};

#endif
