#ifndef FIXED_HPP
#define FIXED_HPP

class FixedEmission : public EmissionFunction {
public:
  FixedEmission(int stateID, int slotID) : EmissionFunction(stateID, slotID), _log_prob(0.0) {}

  virtual bool validParams(Params const & params) const {
    return params.length() == 1 && params[0] >= 0 && params[0] <= 1;
  }

  virtual Params * getParams() const {
    double prob = exp(_log_prob);
    Params * params = new Params(1, &prob);
    params->setFixed(0, true);
    return params;
  }
  
  virtual void setParams(Params const & params) {
    double prob = params[0];
    _log_prob = log(prob);
  }

  virtual double log_probability(Iter const & iter) const {
    return _log_prob;
  }

private:
  double _log_prob;
};

#endif
