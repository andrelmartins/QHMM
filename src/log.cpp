#include "log.hpp"
#include <cstdarg>
#include <cstdio>

#ifdef RQHMM
#define R_USE_C99_IN_CXX
#include <R.h>
#include <R_ext/Print.h>
#endif

/* print error message; saves typing in mains */
void log_msg(const char *warnfmt, ...) {
  va_list args;
 
  va_start(args, warnfmt);

#ifdef RQHMM
  REvprintf(warnfmt, args);
#else
  vfprintf(stderr, warnfmt, args);
#endif
  va_end(args);
}

void log_state_msg(int state, const char *warnfmt, ...) {
  va_list args;

  va_start(args, warnfmt);

#ifdef RQHMM
  REprintf("[%d]: ", state + 1);
  REvprintf(warnfmt, args);
#else
  fprintf(stderr, "[%d]: ", state);
  vfprintf(stderr, warnfmt, args);
#endif
  va_end(args);
}

void log_state_slot_msg(int state, int slot, const char *warnfmt, ...) {
  va_list args;

  va_start(args, warnfmt);

#ifdef RQHMM
  REprintf("[%d][%d]: ", state + 1, slot + 1);
  REvprintf(warnfmt, args);
#else
  fprintf(stderr, "[%d][%d]: ", state, slot);
  vfprintf(stderr, warnfmt, args);
#endif
  va_end(args);
}
