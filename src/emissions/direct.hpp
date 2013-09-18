#ifndef DIRECT_HPP
#define DIRECT_HPP

#include "../base_classes.hpp"
#include "../em_base.hpp"

class DirectEmission : public EmissionFunction {
public:
  DirectEmission(int stateID, int slotID) : EmissionFunction(stateID, slotID), _flip(0.0), _is_log(false) {}

  virtual bool getOption(const char * name, double * out_value) {
    if (!strcmp(name, "flip")) {
      *out_value = _flip;
      return true;
    } else if (!strcmp(name, "is_log")) {
      *out_value = (_is_log ? 1 : 0);
      return true;
    }
    return false;
  }

  virtual bool setOption(const char * name, double value) {
    if (value != 1.0 && value != 0.0) {
      log_msg("invalid option value: %g : should be 0 or 1\n",
	      value);
      return false;
    }

    if (!strcmp(name, "flip")) {
      _flip = value;
      return true;
    } else if (!strcmp(name, "is_log")) {
      _is_log = (value == 1.0);
      return true;
    }
    return false;
  }

  virtual double log_probability(Iter const & iter) const {
    if (_is_log) {
      double log_prob = iter.emission(_slotID);

      if (_flip != 0)
        return log(_flip - exp(log_prob));
      return log_prob;
    } else {
      double prob = iter.emission(_slotID);
      double prob_flipped = _flip * (_flip - prob) + (1.0 - _flip) * prob;

      return log(prob_flipped);
    }
  }

private:
  double _flip;
  bool _is_log;
};

#endif
