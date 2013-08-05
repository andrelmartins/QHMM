#include "rqhmm.hpp"
#include <R_ext/Rdynload.h>
#include <cstdlib>

extern "C" {
  
#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif
  
  // Entry points
  void attribute_hidden rqhmm_register_emission(FuncEntry* emission) {
    static reg_func_t = NULL;
    if (func == NULL)
      func = (reg_func_t) R_GetCCallable("rqhmm", "register_emission");
    return func(emission);
  }
  
  void attribute_hidden rqhmm_register_transition(FuncEntry* transition) {
    static reg_func_t = NULL;
    if (func == NULL)
      func = (reg_func_t) R_GetCCallable("rqhmm", "register_transition");
    return func(transition);
  }

  void attribute_hidden rqhmm_unregister_all(const char * package) {
    static unreg_all_t = NULL;
    if (func == NULL)
      func = (unreg_all_t) R_GetCCallable("rqhmm", "unregister_all");
    return func(package);
  }
}
