#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "func_entry.hpp"

#include <vector>

static std::vector<FuncEntry*> __emissions;
static std::vector<FuncEntry*> __transitions;

extern "C" {
#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attr_hidden __attribute__ ((visibility ("hidden")))
# define attr_default __attribute__ ((visibility ("default")))
#else
# define attr_hidden
# define attr_default
#endif
#define RREGDEF(name)  R_RegisterCCallable("rqhmm", #name, (DL_FUNC) name)

  void attr_hidden register_emission(FuncEntry * function) {
    // TODO: add check to avoid duplicates
    __emissions.push_back(function);
  }
  
  void attr_hidden register_transition(FuncEntry * function) {
    // TODO: add check to avoid duplicates
    __transitions.push_back(function);
  }
  
  void attr_hidden unregister_all(const char * package) {
    std::vector<FuncEntry*>::iterator it;
    
    // drop emissions from package
    it = __emissions.begin();
    while (it != __emissions.end()) {
      FuncEntry * entry = *it;
      if (!strcmp(entry->package, package))
        it = __emissions.erase(it);
      else
        ++it;
    }
    
    // drop transitions from package
    it = __transitions.begin();
    while (it != __emissions.end()) {
      FuncEntry * entry = *it;
      if (!strcmp(entry->package, package))
        it = __emissions.erase(it);
      else
        ++it;
    }
  }
  
  void attr_hidden empty_tables() {
    std::vector<FuncEntry*>::iterator it;
    
    for (it = __emissions.begin(); it != __emissions.end(); it = __emissions.erase(it));
    for (it = __transitions.begin(); it != __transitions.end(); it = __transitions.erase(it));
  }
  
  SEXP attr_hidden fentry_exists(SEXP name, std::vector<FuncEntry*> & table) {
    SEXP result;
    const char * name_ptr;
    
    PROTECT(name = AS_CHARACTER(name));
    PROTECT(result = NEW_LOGICAL(1));
    LOGICAL(result)[0] = FALSE;
    
    name_ptr = CHAR(STRING_ELT(name, 0));
    
    std::vector<FuncEntry*>::iterator it;
    for (it = table.begin(); it != table.end(); ++it)
      if (!strcmp((*it)->name, name_ptr)) {
        LOGICAL(result)[0] = TRUE;
        break;
      }
    
    UNPROTECT(2);
    
    return result;

  }
  
  SEXP rqhmm_transition_exists(SEXP name) {
    return fentry_exists(name, __transitions);
  }
  
  SEXP rqhmm_emission_exists(SEXP name) {
    return fentry_exists(name, __emissions);
  }
  
  // R Entry points
  void attr_default R_init_rqhmm(DllInfo * info) {
    Rprintf("rqhmm init called\n");
    
    // register external routines
    RREGDEF(register_emission);
    RREGDEF(register_transition);
    RREGDEF(unregister_all);
  }
  
  void attr_default R_unload_rqhmm(DllInfo * info) {
    Rprintf("rqhmm unload called\n");
    
    // empty emission & transition tables
    empty_tables();
  }
  
}
