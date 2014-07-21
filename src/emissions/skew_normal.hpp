#ifndef SKEW_NORMAL_HPP
#define SKEW_NORMAL_HPP

#define _USE_MATH_DEFINES
#include <cmath>
#include "../base_classes.hpp"
#include "../em_base.hpp"
#include "../math.hpp"

class SkewNormal : public EmissionFunction {
public:
  SkewNormal(int stateID, int slotID, double location = 0.0, double scale = 1.0, double skew = 1.0) : EmissionFunction(stateID, slotID), _location(location), _scale(scale), _skew(skew), _is_fixed(false), _maxIter(100), _tolerance(1e-6) {}
  
  virtual bool validParams(Params const & params) const {
    return params.length() == 3 && params[1] > 0;
  }
  
  virtual Params * getParams() const {
    double * values = new double[3];
    
    values[0] = _location;
    values[1] = _scale;
    values[2] = _skew;
    
    Params * result = new Params(3, values);
    if (_is_fixed) {
      result->setFixed(0, _is_fixed);
      result->setFixed(1, _is_fixed);
      result->setFixed(2, _is_fixed);
    }
    
    delete[] values;
    
    return result;
  }

  virtual void setParams(Params const & params) {
    _location = params[0];
    _scale = params[1];
    _skew = params[2];
    _is_fixed = params.isFixed(0) || params.isFixed(1) || params.isFixed(2);
    
    update_consts();
  }

  virtual bool getOption(const char * name, double * out_value) {
    if(!strcmp(name, "maxIter")) {
      *out_value = _maxIter;
      return true;
    } else if (!strcmp(name, "tolerance")) {
      *out_value = _tolerance;
      return true;
    }
  }
  
  virtual bool setOption(const char * name, double value) {
    if (!strcmp(name, "maxIter")) {
      int maxIter = (int) value;
      if (maxIter <= 0) {
        log_msg("maxIter must be > 0: %d\n", maxIter);
        return false;
      }
      _maxIter = maxIter;
      return true;
    } else if (!strcmp(name, "tolerance")) {
      if (value < 0) {
        log_msg("tolerance must be >= 0: %g\n", value);
        return false;
      }
      _tolerance = value;
      return true;
    }

    return false;
  }
  
  virtual double log_probability(Iter const & iter) const {
    double x = iter.emission(_slotID);
    return log_pdf(x) + log_skewed_2_cdf(x);
  }

  virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {
    if (_is_fixed)
      return;

    /* auxiliary data */
    struct opt_data data;
    data.sequences = sequences;
    data.group = group;

    double xin[3];
    double fout;
    int fncount = 0;
    int fail = 0;

    /* bounds for L-BFGS-B */
    double lower[3];
    double upper[3];
    int nbd[3];
    int lmm = 5;
    int trace = 0;
    int report = 10;
    double factr = 1e7; // TODO: check if this should be a transformation of the "_tolerance" value
    double pgtol = 0;
    int grcount = 0;
    char msg[60];

    /* initial point */
    xin[0] = _location;
    xin[1] = _scale;
    xin[2] = _skew;

    /* NOTE:
     *
     * Although both Nelder-Mead and L-BFGS-B methods seem to converge
     * equally well on the tests I conducted, L-BFGS-B produces less 
     * function evaluations even when taking into account the additional
     * work to compute the gradients. Moreover, L-BFGS-B presents a natural
     * way to impose the lower bound on scale.
     */
    
    /* optimization via L-BFGS-B */
    lower[0] = -std::numeric_limits<double>::infinity();
    lower[1] = 0.01; /* cannot be zero */
    lower[2] = -std::numeric_limits<double>::infinity();
    upper[0] = std::numeric_limits<double>::infinity();
    upper[1] = std::numeric_limits<double>::infinity();
    upper[2] = std::numeric_limits<double>::infinity();
    nbd[0] = 0; /* unbounded */
    nbd[1] = 1; /* only lower bound */
    nbd[2] = 0; /* unbounded */
    
    QHMM_lbfgsb(3, lmm, xin, lower, upper,
         nbd, &fout, skewnorm_optim, skewnorm_gr, &fail, &data,
         factr, pgtol,
         &fncount, &grcount,
         _maxIter, msg, trace, report);

    if (fail == 0) {
      /* success, propagate values */
      _location = xin[0];
      _scale = xin[1];
      _skew = xin[2];

      update_consts();

      // propagate to other elements in the group
      std::vector<EmissionFunction*>::iterator ef_it;

      for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
        SkewNormal * ef = (SkewNormal*) (*ef_it)->inner();
      
        if (ef != this) {
          ef->_location = _location;
          ef->_scale = _scale;
          ef->_skew = _skew;
          ef->update_consts();
        }
      }
    } else {
      log_state_slot_msg(_stateID, _slotID, "update failed: location = [%g => %g] scale = [%g => %g] skewness = [%g => %g] fail = %d, fout = %g, n.calls = %d\n", 
                         _location, xin[0], 
                         _scale, xin[1],
                         _skew, xin[2],
                         fail, fout, fncount);
    }
  }
  
 private:
  double _location;
  double _scale;
  double _skew;
  bool _is_fixed;
  int _maxIter;
  double _tolerance;

  double _A;  // -log(sqrt(2 pi scale))
  double _B;  // sqrt(2 scale)

  void update_consts() {
    _A = -log(sqrt(2.0 * M_PI * _scale));
    _B = sqrt(2.0 * _scale);
  }

  double log_pdf(double x) const {
    double diff = x - _location;
    
    return _A - (diff * diff) / (2 * _scale);
  }

  double log_skewed_2_cdf(double x) const {
    return log1p(erf(_skew * (x - _location)/_B));
  }

  /* numerical optimization */
  struct opt_data {
    EMSequences * sequences;
    std::vector<EmissionFunction*> * group;

    double f;
  };

  static double full_log_prob(double x, double * params) {
    double location = params[0];
    double scale = params[1];
    double skew = params[3];

    double diff = x - location;
    double A = -log(sqrt(2.0 * M_PI * scale));
    double B = sqrt(2.0 * scale);

    double lpdf = A - (diff * diff) / (2 * scale);
    double ls2cdf = log1p(erf(skew * (x - location) / B));

    return lpdf + ls2cdf;
  }

  static double skewnorm_optim(int n, double *par, void *ex) {
    struct opt_data * data = (struct opt_data*) ex;
    double result = 0;
    
    std::vector<EmissionFunction*>::iterator ef_it;
    for (ef_it = data->group->begin(); ef_it != data->group->end(); ++ef_it) {
      SkewNormal * ef = (SkewNormal*) (*ef_it)->inner();
      PosteriorIterator * post_it = data->sequences->iterator(ef->_stateID, ef->_slotID);
        
      do {
        const double * post_j = post_it->posterior();
        Iter & iter = post_it->iter();
        iter.resetFirst();
        
        for (int j = 0; j < iter.length(); iter.next(), ++j) {
          double x = iter.emission(ef->_slotID);
          
          if (post_j[j] > 0) {
            double logPxi = full_log_prob(x, par);
            result += post_j[j] * logPxi;
          }
        }
      } while (post_it->next());
    }
    
    data->f = -result; /* save for use in skewnorm_gr
                          NOTE: Dependency on this being executed before the
                          gradient */
    return -result;
  }

  /* gr = result vector of length 'n' */
  static void skewnorm_gr(int n, double *par, double *gr, void *ex) {
    struct opt_data * data = (struct opt_data*) ex;
    double result = 0;
    
    int i;
    double f_center = data->f; /* NOTE: Dependency on skewnorm_optim
                                  being executed before the gradient */
    double deriv_epsilon = 1e-6; /* lifted from Phast */
    
    for (i = 0; i < n; ++i) {
      double p_old = par[i];
      double f_step;
      
      /* Forward method (this avoids needing to check lower bounds) */
      par[i] += deriv_epsilon;
      f_step = skewnorm_optim(n, par, ex);
      
      /* save gradient */
      gr[i] = (f_step - f_center) / deriv_epsilon;
      
      /* restore parameter */
      par[i] = p_old;
    }
  }
};

  
#endif
