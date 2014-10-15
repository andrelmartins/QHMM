#ifndef POISSON_H
#define POISSON_H

#define _USE_MATH_DEFINES
#include <cmath>
#include "../base_classes.hpp"
#include "../em_base.hpp"

// Code adapted from http://www.johndcook.com/csharp_log_factorial.html
class LogFactorial {
  
  public:
  
    static double logFactorial(int n) {
      assert(n >= 0);
      
      if (n > 254) {
        int x = n + 1;
        
        return (x - 0.5)*log(x) - x + 0.5*log(2*M_PI) + 1.0/(12.0*x);
      } else
        return logFactorialTable[n];
    }
  
  private:
    static const double logFactorialTable[256];
};

const double LogFactorial::logFactorialTable[256] = { 0.000000000000000,
                0.000000000000000,
                0.693147180559945,
                1.791759469228055,
                3.178053830347946,
                4.787491742782046,
                6.579251212010101,
                8.525161361065415,
                10.604602902745251,
                12.801827480081469,
                15.104412573075516,
                17.502307845873887,
                19.987214495661885,
                22.552163853123421,
                25.191221182738683,
                27.899271383840894,
                30.671860106080675,
                33.505073450136891,
                36.395445208033053,
                39.339884187199495,
                42.335616460753485,
                45.380138898476908,
                48.471181351835227,
                51.606675567764377,
                54.784729398112319,
                58.003605222980518,
                61.261701761002001,
                64.557538627006323,
                67.889743137181526,
                71.257038967168000,
                74.658236348830158,
                78.092223553315307,
                81.557959456115029,
                85.054467017581516,
                88.580827542197682,
                92.136175603687079,
                95.719694542143202,
                99.330612454787428,
                102.968198614513810,
                106.631760260643450,
                110.320639714757390,
                114.034211781461690,
                117.771881399745060,
                121.533081515438640,
                125.317271149356880,
                129.123933639127240,
                132.952575035616290,
                136.802722637326350,
                140.673923648234250,
                144.565743946344900,
                148.477766951773020,
                152.409592584497350,
                156.360836303078800,
                160.331128216630930,
                164.320112263195170,
                168.327445448427650,
                172.352797139162820,
                176.395848406997370,
                180.456291417543780,
                184.533828861449510,
                188.628173423671600,
                192.739047287844900,
                196.866181672889980,
                201.009316399281570,
                205.168199482641200,
                209.342586752536820,
                213.532241494563270,
                217.736934113954250,
                221.956441819130360,
                226.190548323727570,
                230.439043565776930,
                234.701723442818260,
                238.978389561834350,
                243.268849002982730,
                247.572914096186910,
                251.890402209723190,
                256.221135550009480,
                260.564940971863220,
                264.921649798552780,
                269.291097651019810,
                273.673124285693690,
                278.067573440366120,
                282.474292687630400,
                286.893133295426990,
                291.323950094270290,
                295.766601350760600,
                300.220948647014100,
                304.686856765668720,
                309.164193580146900,
                313.652829949878990,
                318.152639620209300,
                322.663499126726210,
                327.185287703775200,
                331.717887196928470,
                336.261181979198450,
                340.815058870798960,
                345.379407062266860,
                349.954118040770250,
                354.539085519440790,
                359.134205369575340,
                363.739375555563470,
                368.354496072404690,
                372.979468885689020,
                377.614197873918670,
                382.258588773060010,
                386.912549123217560,
                391.575988217329610,
                396.248817051791490,
                400.930948278915760,
                405.622296161144900,
                410.322776526937280,
                415.032306728249580,
                419.750805599544780,
                424.478193418257090,
                429.214391866651570,
                433.959323995014870,
                438.712914186121170,
                443.475088120918940,
                448.245772745384610,
                453.024896238496130,
                457.812387981278110,
                462.608178526874890,
                467.412199571608080,
                472.224383926980520,
                477.044665492585580,
                481.872979229887900,
                486.709261136839360,
                491.553448223298010,
                496.405478487217580,
                501.265290891579240,
                506.132825342034830,
                511.008022665236070,
                515.890824587822520,
                520.781173716044240,
                525.679013515995050,
                530.584288294433580,
                535.496943180169520,
                540.416924105997740,
                545.344177791154950,
                550.278651724285620,
                555.220294146894960,
                560.169054037273100,
                565.124881094874350,
                570.087725725134190,
                575.057539024710200,
                580.034272767130800,
                585.017879388839220,
                590.008311975617860,
                595.005524249382010,
                600.009470555327430,
                605.020105849423770,
                610.037385686238740,
                615.061266207084940,
                620.091704128477430,
                625.128656730891070,
                630.172081847810200,
                635.221937855059760,
                640.278183660408100,
                645.340778693435030,
                650.409682895655240,
                655.484856710889060,
                660.566261075873510,
                665.653857411105950,
                670.747607611912710,
                675.847474039736880,
                680.953419513637530,
                686.065407301994010,
                691.183401114410800,
                696.307365093814040,
                701.437263808737160,
                706.573062245787470,
                711.714725802289990,
                716.862220279103440,
                722.015511873601330,
                727.174567172815840,
                732.339353146739310,
                737.509837141777440,
                742.685986874351220,
                747.867770424643370,
                753.055156230484160,
                758.248113081374300,
                763.446610112640200,
                768.650616799717000,
                773.860102952558460,
                779.075038710167410,
                784.295394535245690,
                789.521141208958970,
                794.752249825813460,
                799.988691788643450,
                805.230438803703120,
                810.477462875863580,
                815.729736303910160,
                820.987231675937890,
                826.249921864842800,
                831.517780023906310,
                836.790779582469900,
                842.068894241700490,
                847.352097970438420,
                852.640365001133090,
                857.933669825857460,
                863.231987192405430,
                868.535292100464630,
                873.843559797865740,
                879.156765776907600,
                884.474885770751830,
                889.797895749890240,
                895.125771918679900,
                900.458490711945270,
                905.796028791646340,
                911.138363043611210,
                916.485470574328820,
                921.837328707804890,
                927.193914982476710,
                932.555207148186240,
                937.921183163208070,
                943.291821191335660,
                948.667099599019820,
                954.046996952560450,
                959.431492015349480,
                964.820563745165940,
                970.214191291518320,
                975.612353993036210,
                981.015031374908400,
                986.422203146368590,
                991.833849198223450,
                997.249949600427840,
                1002.670484599700300,
                1008.095434617181700,
                1013.524780246136200,
                1018.958502249690200,
                1024.396581558613400,
                1029.838999269135500,
                1035.285736640801600,
                1040.736775094367400,
                1046.192096209724900,
                1051.651681723869200,
                1057.115513528895000,
                1062.583573670030100,
                1068.055844343701400,
                1073.532307895632800,
                1079.012946818975000,
                1084.497743752465600,
                1089.986681478622400,
                1095.479742921962700,
                1100.976911147256000,
                1106.478169357800900,
                1111.983500893733000,
                1117.492889230361000,
                1123.006317976526100,
                1128.523770872990800,
                1134.045231790853000,
                1139.570684729984800,
                1145.100113817496100,
                1150.633503306223700,
                1156.170837573242400,
            };

class Poisson : public EmissionFunction {
  public:
  Poisson(int stateID, int slotID, double lambda = 1.0) : EmissionFunction(stateID, slotID), _lambda(lambda), _log_lambda(log(lambda)), _is_fixed(false) {}

    virtual bool validParams(Params const & params) const {
      return params.length() == 1 && params[0] > 0;
    }

    virtual Params * getParams() const {
      Params * params = new Params(1, &_lambda);
      params->setFixed(0, _is_fixed);
      return params;
    }

    virtual void setParams(Params const & params) {
      _lambda = params[0];
      _log_lambda = log(_lambda);
      _is_fixed = params.isFixed(0);
    }
  
    virtual double log_probability(Iter const & iter) const {
      int x = (int) iter.emission(_slotID); // cast to integer
      
      // log prob(x) = log( lambda^x exp(-lambda) / x!)
      //             = x log(lambda) - lambda - log(x!)
      if (x == 0)
        return -_lambda - LogFactorial::logFactorial(x);
      else
        return x * _log_lambda - _lambda - LogFactorial::logFactorial(x);
    }
  
    virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {
      if (_is_fixed)
        return;

      // sufficient statistics are the sum of the state posteriors and the sum of the posterior times
      // the observations
      double sum_Pzi = 0;
      double sum_Pzi_xi = 0;
      
      std::vector<EmissionFunction*>::iterator ef_it;
      
      for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
        Poisson * ef = (Poisson*) (*ef_it)->inner();
        PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
        
        do {
          const double * post_j = post_it->posterior();
          Iter & iter = post_it->iter();
          iter.resetFirst();
          
          for (int j = 0; j < iter.length(); iter.next(), ++j) {
            int x = (int) iter.emission(ef->_slotID);
            
            sum_Pzi += post_j[j];
            sum_Pzi_xi += post_j[j] * x;
          }
        } while (post_it->next());
        
        delete post_it;
      }
      
      // use expected counts to estimate parameter value
      _lambda = sum_Pzi_xi / sum_Pzi;
      _log_lambda = log(_lambda);
      
      // propagate to other elements in the group
      for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
        Poisson * ef = (Poisson*) (*ef_it)->inner();
        
        if (ef != this) {
          ef->_lambda = _lambda;
          ef->_log_lambda = _log_lambda;
        }
      }
    }
  
  private:
    double _lambda;
    double _log_lambda;
    bool _is_fixed;
};

class PoissonCovar : public EmissionFunction {
  public:
  PoissonCovar(int stateID, int slotID, int covar_slot = 0) : EmissionFunction(stateID, slotID), _covar_slot(covar_slot) {}
    
    virtual double log_probability(Iter const & iter) const {
      int x = (int) iter.emission(_slotID); // cast to integer
      double lambda = iter.covar(_covar_slot);
      
      // log prob(x) = log( lambda^x exp(-lambda) / x!)
      //             = x log(lambda) - lambda - log(x!)
      if (x == 0)
        return - lambda - LogFactorial::logFactorial(x);
      else
        return x * log(lambda) - lambda - LogFactorial::logFactorial(x);
    }
  
    virtual bool setCovarSlots(int * slots, int length) {
      if (length != 1)
        return false;
      if (*slots < 0)
        return false;
      _covar_slot = *slots;
      return true;
    }
  
  private:
    int _covar_slot;
};

class PoissonScaledCovar : public EmissionFunction {
public:
  
  PoissonScaledCovar(int stateID, int slotID, int covar_slot = 0, double scale = 1) : EmissionFunction(stateID, slotID), _covar_slot(covar_slot), _scale(scale) {
    _is_fixed = false;
    _lower_bound = 0;
    _upper_bound = std::numeric_limits<double>::infinity();
    _pseudo_num = 0;
    _pseudo_denom = 0;
  }
  
  virtual bool validParams(Params const & params) const {
    return params.length() == 1 && params[0] > 0;
  }
  
  virtual Params * getParams() const {
    Params * params = new Params(1, &_scale);
    params->setFixed(0, _is_fixed);
    return params;
  }
  
  virtual void setParams(Params const & params) {
    _scale = params[0];
    _is_fixed = params.isFixed(0);
  }
  
  virtual bool getOption(const char * name, double * out_value) {
    if (!strcmp(name, "lower_bound")) {
      *out_value = _lower_bound;
      return true;
    } else if (!strcmp(name, "upper_bound")) {
      *out_value = _upper_bound;
      return true;
    } else if (!strcmp(name, "pseudo_num")) {
      *out_value = _pseudo_num;
      return true;
    } else if (!strcmp(name, "pseudo_denom")) {
      *out_value = _pseudo_denom;
      return true;
    }
    return false;
  }
  
  virtual bool setOption(const char * name, double value) {
    if (!strcmp(name, "lower_bound")) {
      if (value < 0 || value >= _upper_bound) {
        log_msg("invalid lower_bound: %g : must be between > 0 and <= %g\n",
                value, _upper_bound);
        return false;
      }
      _lower_bound = value;
      return true;
    } else if (!strcmp(name, "upper_bound")) {
      if (value <= _lower_bound) {
        log_msg("invalid upper_bound: %g : must be > %g\n",
                _lower_bound);
        return false;
      }
      _upper_bound = value;
    } else if (!strcmp(name, "pseudo_num")) {
      if (value < 0) {
        log_msg("invalid pseudo_num: %g : must be >= 0\n",
                value);
        return false;
      }
      _pseudo_num = value;
      return true;
    } else if (!strcmp(name, "pseudo_denom")) {
      if (value < 0) {
        log_msg("invalid pseudo_denom: %g : must be >= 0\n",
                value);
        return false;
      }
      _pseudo_denom = value;
      return true;
    }
    return false;
  }

  virtual double log_probability(Iter const & iter) const {
    int x = (int) iter.emission(_slotID); // cast to integer
    double lambda = iter.covar(_covar_slot) * _scale;
    
    // log prob(x) = log( lambda^x exp(-lambda) / x!)
    //             = x log(lambda) - lambda - log(x!)
    if (x == 0)
      return - lambda - LogFactorial::logFactorial(x);
    else
      return x * log(lambda) - lambda - LogFactorial::logFactorial(x);
  }
  
  virtual bool setCovarSlots(int * slots, int length) {
    if (length != 1)
      return false;
    if (*slots < 0)
      return false;
    _covar_slot = *slots;
    return true;
  }

  virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {
    if (_is_fixed)
      return;

    // sufficient statistics are the sum of the state posteriors times the lambdas and the sum of the posterior times
    // the observations
    double sum_Pzi_lambda_i = 0;
    double sum_Pzi_xi = 0;
      
    std::vector<EmissionFunction*>::iterator ef_it;
    
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      PoissonScaledCovar * ef = (PoissonScaledCovar*) (*ef_it)->inner();
      PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
      
      do {
        const double * post_j = post_it->posterior();
        Iter & iter = post_it->iter();
        iter.resetFirst();
        
        for (int j = 0; j < iter.length(); iter.next(), ++j) {
          int x = (int) iter.emission(ef->_slotID);
          double lambda_i = iter.covar(_covar_slot);
          
          sum_Pzi_lambda_i += post_j[j] * lambda_i;
          sum_Pzi_xi += post_j[j] * x;
        }
      } while (post_it->next());
      
      delete post_it;
    }
    
    // update parameter
    double scale = (_pseudo_num + sum_Pzi_xi) / (_pseudo_denom + sum_Pzi_lambda_i);

    if (scale > _lower_bound && scale < _upper_bound) {
      _scale = scale;

      // propagate to other elements in the group
      for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
        PoissonScaledCovar * ef = (PoissonScaledCovar*) (*ef_it)->inner();
        
        if (ef != this)
          ef->_scale = _scale;
      }

    } else
      log_state_slot_msg(_stateID, _slotID, "scale update passed bounds: %g [%g, %g]\n", scale, _lower_bound, _upper_bound);
  }
  
private:
  int _covar_slot;
  double _scale;
  double _lower_bound;
  double _upper_bound;
  double _pseudo_num;
  double _pseudo_denom;
  bool _is_fixed;
};


class PoissonScaled : public EmissionFunction {
  public:
    PoissonScaled(int stateID, int slotID, double lambda = 1.0) : EmissionFunction(stateID, slotID), _lambda(lambda), _scale(1.0), _scale_lambda(lambda), _log_scale_lambda(log(lambda)), _is_fixed(false) {}

    virtual bool validParams(Params const & params) const {
      return params.length() == 1 && params[0] > 0;
    }

    virtual Params * getParams() const {
      Params * params = new Params(1, &_lambda);
      params->setFixed(0, _is_fixed);
      return params;
    }

    virtual void setParams(Params const & params) {
      update_lambda(params[0]);
      _is_fixed = params.isFixed(0);
    }

    virtual bool getOption(const char * name, double * out_value) {
      if (!strcmp(name, "scale")) {
        *out_value = (double) _scale;
        return true;
      }
      return false;
    }

    virtual bool setOption(const char * name, double value) {
      if (!strcmp(name, "scale")) {
        if (value <= 0) {
          log_msg("scale must be > 0: %g\n", value);
          return false;
        }
        update_scale(value);
        return true;
      }
      return false;
    }
    
    virtual double log_probability(Iter const & iter) const {
      int x = (int) iter.emission(_slotID); // cast to integer
      
      // log prob(x) = log( (scale*lambda)^x exp(-scale*lambda) / x!)
      //             = x log(scale*lambda) - scale*lambda - log(x!)
      if (x == 0)
        return -_scale_lambda - LogFactorial::logFactorial(x);
      else
        return x * _log_scale_lambda - _scale_lambda - LogFactorial::logFactorial(x);
    }

    virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {
      if (_is_fixed)
        return;

      // sufficient statistics are the sum of 'scaled' the state posteriors and the sum of the posterior times
      // the observations
      double sum_scale_Pzi = 0;
      double sum_Pzi_xi = 0;
      
      std::vector<EmissionFunction*>::iterator ef_it;
      
      for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
        PoissonScaled * ef = (PoissonScaled*) (*ef_it)->inner();
        PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
        double scale_i = ef->_scale;
        
        do {
          const double * post_j = post_it->posterior();
          Iter & iter = post_it->iter();
          iter.resetFirst();
          
          for (int j = 0; j < iter.length(); iter.next(), ++j) {
            int x = (int) iter.emission(ef->_slotID);
            
            sum_scale_Pzi += post_j[j] * scale_i;
            sum_Pzi_xi += post_j[j] * x;
          }
        } while (post_it->next());
        
        delete post_it;
      }
      
      // if no data was found, don't update the parameters
      if (sum_scale_Pzi == 0) {
        log_state_slot_msg(_stateID, _slotID, "no data, skipping update\n");
        return;
      }
      
      // use expected counts to estimate parameter value
      update_lambda(sum_Pzi_xi / sum_scale_Pzi);
      
      // propagate to other elements in the group
      for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
        PoissonScaled * ef = (PoissonScaled*) (*ef_it)->inner();
        
        if (ef != this)
          ef->update_lambda(_lambda);
      }
    }

  private:
    double _lambda;
    double _scale;
    double _scale_lambda;
    double _log_scale_lambda;
    bool _is_fixed;
        
    void update_scale_lambda(double lambda, double scale) {
      _scale = scale;
      _lambda = lambda;
      _scale_lambda = _scale * _lambda;
      _log_scale_lambda = log(scale) + log(lambda);
    }
    
    void update_scale(double scale) {
      update_scale_lambda(scale, _lambda);
    }
    
    void update_lambda(double lambda) {
      update_scale_lambda(_scale, lambda);
    }
};

#endif
