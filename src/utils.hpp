#ifndef UTILS_HPP
#define UTILS_HPP

#include "fn_name.hpp"

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"

#include <memory>
#include <string>
#include <vector>

std::vector<size_t>
get_random_indices(size_t total_size, size_t count);

void
wait_till_parent_ends();

void
bind_to_parent();

void
make_pipe(const std::string& pipepath);

std::vector<std::string>
read_pipe_lines(const std::string& pipepath, const std::string& end_line);

std::string
read_pipe(const std::string& pipepath);

void
confirm_pipe(const std::string& pipepath);

void
fill_bfs(const char* seq,
         size_t seq_len,
         unsigned hash_num,
         const std::vector<unsigned>& k_values,
         unsigned kmer_threshold,
         std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
         std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs);

inline void
fill_bfs(const std::string& seq,
         unsigned hash_num,
         const std::vector<unsigned>& k_values,
         unsigned kmer_threshold,
         std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
         std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs)
{
  fill_bfs(
    seq.c_str(), seq.size(), hash_num, k_values, kmer_threshold, cbfs, bfs);
}

#endif