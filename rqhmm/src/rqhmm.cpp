#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" {
#define RREGDEF(name)  R_RegisterCCallable("rqhmm", #name, (DL_FUNC) name)
  
  // R Entry points
  void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
  __attribute__ ((visibility ("default")))
#endif
  R_init_rqhmm(DllInfo * info) {
    Rprintf("rqhmm init called\n");
    
    // register external routines
    
  }

  void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
  __attribute__ ((visibility ("default")))
#endif
  R_unload_rqhmm(DllInfo * info) {
    Rprintf("rqhmm unload called\n");
    
  }
  
  SEXP rqhmm_transition_exists(SEXP name) {
    SEXP result;
    
    PROTECT(result = NEW_LOGICAL(1));
    
    LOGICAL(result)[0] = FALSE;
    
    UNPROTECT(1);
    
    return result;
  }
  
  SEXP rqhmm_emission_exists(SEXP name) {
    SEXP result;
    
    PROTECT(result = NEW_LOGICAL(1));
    
    LOGICAL(result)[0] = FALSE;
    
    UNPROTECT(1);
    
    return result;
  }
  
  
}
