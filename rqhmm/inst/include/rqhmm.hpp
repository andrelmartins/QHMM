#ifndef RQHMM_HPP
#define RQHMM_HPP

#include "func_entry.hpp"

extern "C" {

  void rqhmm_register_emission(FuncEntry* emission);
  void rqhmm_register_transition(FuncEntry* transition);
  void rqhmm_unregister_all(const char * package);
}

#endif
