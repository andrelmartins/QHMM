#include "utils.hpp"
#include <algorithm>

static bool is_block_state(const int state, const std::vector<int> * block_states) {
  std::vector<int>::const_iterator res = std::find(block_states->begin(), block_states->end(), state);
  return res != block_states->end();
}

std::vector<block_t> * path_blocks(const std::vector<int> * path, const std::vector<int> * block_states) {
  std::vector<block_t> * result = new std::vector<block_t>();
  bool in_block = false;
  int block_start;
  std::vector<int>::const_iterator it;
  
  for (it = path->begin(); it != path->end(); ++it) {
    bool block_state = is_block_state(*it, block_states);
    
    if (in_block && !block_state) {
      block_t block;
      block.start = block_start;
      block.end = it - path->begin() - 1;
      result->push_back(block);
      in_block = false;
    } else if (!in_block && block_state) {
      in_block = true;
      block_start = it - path->begin();
    }
  }
  
  // last block
  if (in_block) {
    block_t block;
    block.start = block_start;
    block.end = path->size() - 1;
    result->push_back(block);
  }
  
  return result;
}

std::vector<block_t> * path_blocks(const std::vector<int> * path, const std::vector<int> * start_states, const std::vector<int> * middle_states, const std::vector<int> * end_states) {
  std::vector<block_t> * result = new std::vector<block_t>();
  bool in_block = false;
  bool seen_mid = false;
  bool seen_end = false;
  int block_start;
  std::vector<int>::const_iterator it;
  
  for (it = path->begin(); it != path->end(); ++it) {
    if (in_block) {
      if (is_block_state(*it, middle_states)) {
        seen_mid = true;
        continue;
      }
      
      if ((seen_mid || seen_end) && is_block_state(*it, end_states)) {
	seen_end = true;
	continue;
      }

      if (seen_end) {
        // good block
        block_t block;
        block.start = block_start;
        block.end = it - path->begin() - 1;
        result->push_back(block);
        in_block = false;
        seen_mid = false;
        seen_end = false;
        // let things go through so that we can start 
        // matching a new start
      }
      
      // restart match
      if (is_block_state(*it, start_states)) {
        block_start = it - path->begin();
        in_block = true;
        seen_mid = false;
        continue;
      }
      
      // mismatch
      in_block = false;
      seen_mid = false;
    } else if (is_block_state(*it, start_states)) {
      block_start = it - path->begin();
      in_block = true;
    }
  }
  
  // no need to terminate last block since it's a mismatch (no end state)!
  
  return result;
}

std::vector<block_t> * path_blocks_seq(const std::vector<int> * path, const std::vector<int> * block_states) {
  std::vector<block_t> * result = new std::vector<block_t>();
  bool in_block = false;
  int block_start;
  unsigned int block_index = 0;
  std::vector<int>::const_iterator it;
  
  for (it = path->begin(); it != path->end(); ++it) {
    if (in_block) {
      if (*it == (*block_states)[block_index])
        continue;
      // end ?
      if (block_index == block_states->size() - 1) {
        block_t block;
        block.start = block_start;
        block.end = it - path->begin() - 1;
        result->push_back(block);
        in_block = false;
        block_index = 0;
        continue;
      }
      // try next
      ++block_index;
      if (*it == (*block_states)[block_index])
        continue;
      // mismatch
      in_block = false;
      block_index = 0;
    } else if (*it == (*block_states)[0]) {
      block_start = it - path->begin();
      in_block = true;
    }
  }
  
  // last block
  if (in_block && block_index == block_states->size() - 1) {
    block_t block;
    block.start = block_start;
    block.end = path->size()  - 1;
    result->push_back(block);
  }
  
  return result;
}
