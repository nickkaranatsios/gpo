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

#include "gpo/common.h"
#include "gpo/cache.h"

// If lseek64 is not available then fall back on lseek().
#if GPO_LSEEK64_AVAILABLE == 0
#define lseek64 lseek
#endif

using namespace GPO;

Cache::Cache() {
  backing_filename_ = 0;
  fd_ = 0;
  block_size_ = 0;
  cache_blocks_ = 0;
  cache_blocks_m1_ = 0;
  cache_ = 0;
  tag_ = 0;
  group_size_ = 0;
  access_count_ = 0;
  miss_count_ = 0;
}

Cache::~Cache() {
  Reset();
}

void Cache::Initialize(size_t block_size, size_t cache_blocks,
                        const char *backing_filename) {
  Reset();

  // Sanity check arguments. 'cache_blocks' must be a power of 2.
  CHECK((cache_blocks & (cache_blocks-1)) == 0);

  // Open backing file, make sure it's writable. Assign each opened file a
  // unique number to prevent collisions in case the caller uses the same
  // filename for different caches.
  static int unique_number = 0;
  asprintf(&backing_filename_, "CACHE_%s.%d", backing_filename, unique_number);
  unique_number++;
  fd_ = open(backing_filename_, O_RDWR | O_TRUNC | O_CREAT | O_LARGEFILE, 0666);
  if (fd_ == -1) {
    GPO_PANIC("Can't open backing file %s (%s)", backing_filename_, strerror(errno));
  }

  // Determine the number of blocks to read/write at once (the group size).
  // The group size is > 1 on the assumption that the user will be doing roughly
  // sequential access and close-by blocks will also need to be brought in to
  // the cache at roughly the same time, so it's faster to do one read/write
  // than many small ones.
  group_size_ = cache_blocks / 8;
  if (group_size_ < 1) group_size_ = 1;

  // Allocate cache memory and tag arrays.
  block_size_ = block_size;
  cache_blocks_ = cache_blocks;
  cache_blocks_m1_ = cache_blocks-1;
  cache_ = (char*) malloc(cache_blocks_ * block_size_);
  if (!cache_) {
    GPO_PANIC("Out of memory");
  }
  tag_ = new size_t[cache_blocks_];
  if (!tag_) {
    GPO_PANIC("Out of memory");
  }

  // Load the tag array with sequential block indexes, as we initially load the
  // cache with zero.
  for (size_t i = 0; i < cache_blocks_; i++) {
    tag_[i] = i;
  }

  // Set the cache to zero
  memset(cache_, 0, cache_blocks_ * block_size_);
}

void Cache::SetZero() {
  if (backing_filename_ && tag_ && cache_) {
    CHECK(ftruncate(fd_, 0) == 0);
    for (size_t i = 0; i < cache_blocks_; i++) {
      tag_[i] = i;
    }
    memset(cache_, 0, cache_blocks_ * block_size_);
  }
}

void Cache::Reset() {
  // Deallocate resources if necessary
  if (backing_filename_) {
    close(fd_);
    unlink(backing_filename_);
    free(backing_filename_);
  }
  if (cache_) {
    free(cache_);
  }
  if (tag_) {
    delete[] tag_;
  }

  // Reset object state
  backing_filename_ = 0;
  fd_ = 0;
  block_size_ = 0;
  cache_blocks_ = 0;
  cache_blocks_m1_ = 0;
  cache_ = 0;
  tag_ = 0;
  group_size_ = 0;
  access_count_ = 0;
  miss_count_ = 0;
}

void Cache::LoadBlock(size_t block_number) {
  // If the block is already in memory, nothing to do
  size_t cache_ofs = block_number & cache_blocks_m1_;
  if (tag_[cache_ofs] == block_number) {
    return;
  }

  // Otherwise we write the colliding block group to disk and then load a new
  // block group in its place.
  size_t block_number_lo = block_number & (~(group_size_-1));
  size_t cache_ofs_lo = block_number_lo & cache_blocks_m1_;
  char *cache_address_lo = cache_ + block_size_ * cache_ofs_lo;
  off64_t file_offset_W = (off64_t)block_size_ * (off64_t)(tag_[cache_ofs_lo]);
  off64_t file_offset_R = (off64_t)block_size_ * (off64_t)block_number_lo;
  CHECK(lseek64(fd_, file_offset_W, SEEK_SET) == file_offset_W);
  CHECK(write(fd_, cache_address_lo, group_size_*block_size_) == (ssize_t)(group_size_*block_size_));
  CHECK(lseek64(fd_, file_offset_R, SEEK_SET) == file_offset_R);
  ssize_t read_status = read(fd_, cache_address_lo, group_size_*block_size_);
  if (read_status == 0) {
    // We read off the end of the file, which means that the loaded block
    // should be all zeros
    memset(cache_address_lo, 0, group_size_*block_size_);
  } else {
    // If we get data, make sure it's a whole block. Partial blocks should
    // never be present in the backing file.
    CHECK(read_status == (ssize_t)(group_size_*block_size_));
  }
  for (ssize_t i = 0; i < group_size_; i++) {
    tag_[cache_ofs_lo+i] = block_number_lo + i;
  }
  miss_count_ += group_size_;
}
