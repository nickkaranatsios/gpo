// Copyright 2007 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Author: Russell L. Smith

// A cache to allow access to blocks in very large (e.g. more than 100 Gb) block
// arrays. It is designed to give fast block access when the block is resident
// in memory, so a direct-mapped (non-associative) cache structure is used.
// In other words, a high cache hit rate is assumed, e.g. arising from a
// roughly sequential access pattern.
//
// On 32 bit machines the block index must be 32 bits, so to get addressing
// spaces larger than 4 Gb you must use block sizes greater than 1. Typical
// usage is with a block size of e.g. 100 bytes, which gives a maximimum
// addressable space of 400 Gb.

#ifndef __GPO_CACHE_H__
#define __GPO_CACHE_H__

#include "gpo/common.h"


// Define the following macro to collect Hessian cache hit/miss statistics.
// This will slow the cache down a little bit, especially for hits.
#define __CACHE_COLLECT_STATS__


namespace GPO {

class Cache {
 public:
  Cache();
  ~Cache();

  // Initialize the cache. The cache_blocks must be a power of two. This sets
  // all cached blocks to zero. The backing file will be truncated by this call,
  // and the destructor will delete this file.
  void Initialize(size_t block_size, size_t cache_blocks,
                   const char *backing_filename);

  // Reset all cached blocks to zero.
  void SetZero();

  // Make sure that the given block is mapped into memory somewhere and return
  // a pointer to it. Flush blocks to disk if necessary to make room. Callers
  // must not keep block pointers beyond the next call to Access() as the
  // underlying data may be moved away.
  void *Access(size_t block_number) {
    #ifdef __CACHE_COLLECT_STATS__
      access_count_++;
    #endif
    size_t b = block_number & cache_blocks_m1_;
    //@@@@ size_t b = block_number % cache_blocks_;
    // If the block is not in memory then make it available.
    if (tag_[b] != block_number) {
      LoadBlock(block_number);
    }
    // The block is in memory, return a pointer to it.
    return cache_ + block_size_ * b;
  }

  // Reset the state of the object to what it is just after the constructor
  // is called.
  void Reset();

  // Return the cache size.
  size_t GetCacheSize() const {
    return cache_blocks_;
  }

  // Return the block group size.
  size_t GetGroupSize() const {
    return group_size_;
  }

  // Return the cache miss rate, if statistics collection is enabled.
  double GetMissRate() const {
    return double(miss_count_) / double(access_count_);
  }

 private:
  char *backing_filename_;      // Name of the backing file
  int fd_;                      // File descriptor of backing file
  size_t block_size_;           // Size in byte of each block
  size_t cache_blocks_;         // Size of the cache, in blocks (power of 2)
  size_t cache_blocks_m1_;      // cache_blocks_ - 1
  char *cache_;                 // A big array of cached blocks
  size_t *tag_;                 // Tag array: block numbers for cached blocks
  size_t group_size_;           // Number of blocks to R/W at once
  size_t access_count_;         // Number of times Access() called
  size_t miss_count_;           // Number of times blocks read by LoadBlock()

  // Make sure the given block is loaded into memory.
  void LoadBlock(size_t block_number);
};

}; // namespace GPO

#undef __CACHE_COLLECT_STATS__
#endif
