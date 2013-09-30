#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>

typedef struct {
  int start;
  int end;
} block_t;

/*
 path_blocks: find contiguous sub-sequences of block states in path
              if start, end and middle block sets are given, then a block must start with one of the start states,
              be filled with any sequence of middle states and end with one or more of the end states
 */
std::vector<block_t> * path_blocks(const std::vector<int> * path, const std::vector<int> * block_states);
std::vector<block_t> * path_blocks(const std::vector<int> * path, const std::vector<int> * start_states, const std::vector<int> * middle_states, const std::vector<int> * end_states);

/* path_blocks_seq: find contiguous sub-sequences of blocks such that each of the block states occurs one or more times within the block.
 
  i.e.: if block_states = {a, b, c} then matches are to reg expr. [a]+[b]+[c]+
*/
std::vector<block_t> * path_blocks_seq(const std::vector<int> * path, const std::vector<int> * block_states);

#endif
