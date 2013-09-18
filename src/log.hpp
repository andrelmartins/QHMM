#ifndef LOG_HPP
#define LOG_HPP

void log_msg(const char *warnfmt, ...);
void log_state_msg(int state, const char *warnfmt, ...);
void log_state_slot_msg(int state, int slot, const char *warnfmt, ...);

#endif
