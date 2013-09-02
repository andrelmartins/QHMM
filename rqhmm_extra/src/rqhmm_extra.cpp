#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <rqhmm_stubs.cpp>
#include <rqhmm.hpp>

extern "C" {
#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attr_hidden __attribute__ ((visibility ("hidden")))
# define attr_default __attribute__ ((visibility ("default")))
#else
# define attr_hidden
# define attr_default
#endif

  // R Entry points
  void attr_default R_init_rqhmm_extra(DllInfo * info) {
    Rprintf("rqhmm.extra init called\n");

  }

  void attr_default R_unload_rqhmm_extra(DllInfo * info) {
    Rprintf("rqhmm.extra unload called\n");
    rqhmm_unregister_all("rqhmm_extra");
  }
}
