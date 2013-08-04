#ifndef FUNC_ENTRY_HPP
#define FUNC_ENTRY_HPP

class FuncEntry {
public:
  const char * package;
  const char * name;
};


typedef void (*reg_func_t)(FuncEntry*);
typedef void (*unreg_all_t)(const char*);

#endif
