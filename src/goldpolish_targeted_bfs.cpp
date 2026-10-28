#include "mappings.hpp"
#include "seqindex.hpp"
#include "utils.hpp"

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static const unsigned MX_THRESHOLD_MIN = 1;
static const unsigned MX_THRESHOLD_MAX = 30;
static const std::string BATCH_NAME_INPUT_PIPE = "batch_name_input";
static const std::string BATCH_TARGET_IDS_INPUT_READY_PIPE =
  "batch_target_ids_input_ready";
static const std::string TARGET_IDS_INPUT_PIPE = "target_ids_input";
static const std::string BFS_READY_PIPE = "bfs_ready";
static const std::string SEPARATOR = "-";
static const std::string BF_EXTENSION = ".bf";
static const std::string END_SYMBOL = "x";

int
mappings_bases_to_kmer_threshold(const unsigned long mappings_bases)
{
  static const double a = 4.66943;
  static const double b = 2.11391e-07;
  static const int max_kmer_threshold = 13;
  const int kmer_threshold = int(std::round(a + double(mappings_bases) * b));
  return std::min(kmer_threshold, max_kmer_threshold);
}

void
serve_batch(const SeqIndex& target_seqs_index,
            const SeqIndex& mapped_seqs_index,
            const AllMappings& all_mappings,
            const size_t cbf_bytes,
            const size_t bf_bytes,
            const std::string& batch_name,
            const std::string& target_ids_input_pipe,
            const std::string& bfs_ready_pipe,
            const std::vector<std::string>& bf_names, // NOLINT
            const unsigned hash_num,
            const std::vector<unsigned>& k_values, // NOLINT
            const double subsample_max_mapped_seqs_per_target_10kbp)
{
  // Initialize Bloom filters
  std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>> cbfs;
  std::vector<std::unique_ptr<btllib::KmerBloomFilter>> bfs;
  for (const auto k : k_values) {
    cbfs.push_back(std::unique_ptr<btllib::KmerCountingBloomFilter8>(
      new btllib::KmerCountingBloomFilter8(cbf_bytes, hash_num, k)));
    bfs.push_back(std::unique_ptr<btllib::KmerBloomFilter>(
      new btllib::KmerBloomFilter(bf_bytes, hash_num, k)));
  }

  // Set Bloom filter paths
  std::vector<std::string> bf_full_names;
  for (const auto& bf_name : bf_names) {
    // NOLINTNEXTLINE
    bf_full_names.push_back(batch_name + SEPARATOR + bf_name);
  }

  for (const auto& id_line : read_pipe_lines(target_ids_input_pipe, END_SYMBOL)) {
    std::istringstream id_stream(id_line);
    std::string target_seq_id;
    if (!(id_stream >> target_seq_id)) {
      continue;
    }
    const auto target_seq_len = target_seqs_index.get_seq_len(target_seq_id);

    const auto& mappings = all_mappings.get_mappings(target_seq_id);
    if (mappings.empty()) {
      continue;
    }
    const auto mappings_num = mappings.size();
    const auto mappings_num_max =
      std::remove_const<decltype(mappings_num)>::type(
        double(target_seq_len) * subsample_max_mapped_seqs_per_target_10kbp /
        10'000.0);
    const auto mappings_num_adjusted = std::min(mappings_num, mappings_num_max);

    std::vector<std::tuple<SeqId, size_t>> mappings_phred;
    for (const auto& mapped_id : mappings) {
      const auto mapped_seq_phred = mapped_seqs_index.get_phred_avg(mapped_id);
      mappings_phred.emplace_back(mapped_id, mapped_seq_phred);
    }

    std::sort(mappings_phred.begin(),
              mappings_phred.end(),
              [](const auto& a, const auto& b) {
                const auto& [a_id, a_phred] = a;
                const auto& [b_id, b_phred] = b;
                return (a_phred > b_phred) ||
                       (a_phred == b_phred && a_id < b_id);
              });

    unsigned long mappings_bases = 0;

    for (size_t i = 0; i < mappings_num_adjusted; i++) {
      const auto& [mapped_id, mapped_seq_phred] = mappings_phred[i];
      const auto mapped_seq_len = mapped_seqs_index.get_seq_len(mapped_id);
      mappings_bases += mapped_seq_len;
    }
    const auto kmer_threshold =
      mappings_bases_to_kmer_threshold(mappings_bases);
    btllib::check_error(kmer_threshold <= 0,
                        FN_NAME + ": k-mer threshold must be >0.");

    // NOLINTNEXTLINE(google-readability-braces-around-statements,hicpp-braces-around-statements,readability-braces-around-statements)
    for (size_t i = 0; i < mappings_num_adjusted; i++) {
      const auto& [mapped_id, mapped_seq_phred] = mappings_phred[i];
      const auto [seq, seq_len] = mapped_seqs_index.get_seq<1>(mapped_id);
      fill_bfs(seq, seq_len, hash_num, k_values, kmer_threshold, cbfs, bfs);
    }
  }

  for (size_t i = 0; i < bfs.size(); i++) {
    bfs[i]->save(bf_full_names[i]);
  }

  confirm_pipe(bfs_ready_pipe);

  std::remove(target_ids_input_pipe.c_str());
  std::remove(bfs_ready_pipe.c_str());
}

bool
process_batch_name(const SeqIndex& target_seqs_index,
                   const SeqIndex& mapped_seqs_index,
                   const AllMappings& all_mappings,
                   const size_t cbf_bytes,
                   const size_t bf_bytes,
                   const std::string& batch_name_input_pipe,
                   const std::string& batch_target_ids_input_ready_pipe,
                   const std::string& target_ids_input_pipe,
                   const std::string& bfs_ready_pipe,
                   const unsigned hash_num,
                   const std::vector<unsigned>& k_values,    // NOLINT
                   const std::vector<std::string>& bf_names, // NOLINT
                   const double subsample_max_mapped_seqs_per_target_10kbp)
{
  const auto batch_name = read_pipe(batch_name_input_pipe);
  if (batch_name.empty() || batch_name == END_SYMBOL) {
    return false;
  }

  const auto batch_target_ids_input_pipe =
    batch_name + SEPARATOR + target_ids_input_pipe;
  const auto batch_bfs_ready_pipe = batch_name + SEPARATOR + bfs_ready_pipe;

  make_pipe(batch_target_ids_input_pipe);
  make_pipe(batch_bfs_ready_pipe);

  confirm_pipe(batch_target_ids_input_ready_pipe);

#pragma omp task firstprivate(                                                 \
    batch_name, batch_target_ids_input_pipe, batch_bfs_ready_pipe)             \
  shared(                                                                      \
      target_seqs_index, mapped_seqs_index, all_mappings, bf_names, k_values)
  serve_batch(target_seqs_index,
              mapped_seqs_index,
              all_mappings,
              cbf_bytes,
              bf_bytes,
              batch_name,
              batch_target_ids_input_pipe,
              batch_bfs_ready_pipe,
              bf_names,
              hash_num,
              k_values,
              subsample_max_mapped_seqs_per_target_10kbp);

  return true;
}

void
serve(const SeqIndex& target_seqs_index,
      const SeqIndex& mapped_seqs_index,
      const AllMappings& all_mappings,
      const size_t cbf_bytes,
      const size_t bf_bytes,
      const std::string& batch_name_input_pipe,
      const std::string& batch_target_ids_input_ready_pipe,
      const std::string& target_ids_input_pipe,
      const std::string& bfs_ready_pipe,
      const unsigned hash_num,
      const std::vector<unsigned>& k_values, // NOLINT
      const double subsample_max_mapped_seqs_per_target_10kbp)
{
  make_pipe(batch_name_input_pipe);
  make_pipe(batch_target_ids_input_ready_pipe);

  std::vector<std::string> bf_names;
  for (const auto k : k_values) {
    // NOLINTNEXTLINE
    bf_names.push_back("k" + std::to_string(k) + BF_EXTENSION);
  }

  btllib::log_info(FN_NAME + ": Accepting batch names at " +
                   batch_name_input_pipe);

#pragma omp parallel
#pragma omp single
  while (process_batch_name(target_seqs_index,
                            mapped_seqs_index,
                            all_mappings,
                            cbf_bytes,
                            bf_bytes,
                            batch_name_input_pipe,
                            batch_target_ids_input_ready_pipe,
                            target_ids_input_pipe,
                            bfs_ready_pipe,
                            hash_num,
                            k_values,
                            bf_names,
                            subsample_max_mapped_seqs_per_target_10kbp)) {
  }

  std::remove(batch_name_input_pipe.c_str());
  std::remove(batch_target_ids_input_ready_pipe.c_str());

  btllib::log_info(FN_NAME + ": Targeted BF builder done!");
}

int
main(int argc, char** argv)
{
  // NOLINTNEXTLINE(readability-magic-numbers,cppcoreguidelines-avoid-magic-numbers)
  btllib::check_error(argc < 8, "Wrong args.");

  bind_to_parent();

  std::vector<unsigned> k_values;

  int arg = 1;
  auto* const target_seqs_filepath = argv[arg++];
  auto* const target_seqs_index_filepath = argv[arg++];
  auto* const mappings_filepath = argv[arg++];
  auto* const mapped_seqs_filepath = argv[arg++];
  auto* const mapped_seqs_index_filepath = argv[arg++];
  const auto mx_max_mapped_seqs_per_target_10kbp = std::stod(argv[arg++]);
  const auto subsample_max_mapped_seqs_per_target_10kbp =
    std::stod(argv[arg++]);
  const auto threads = std::stoi(argv[arg++]);
  while (arg < argc) {
    k_values.push_back(std::stoi(argv[arg++]));
  }

  size_t cbf_bytes = 10ULL * 1024ULL * 1024ULL; // NOLINT
  size_t bf_bytes = 512ULL * 1024ULL;           // NOLINT
  unsigned hash_num = 4;                        // NOLINT

  omp_set_nested(1);
  omp_set_num_threads(int(threads));

  SeqIndex target_seqs_index(target_seqs_index_filepath, target_seqs_filepath);
  SeqIndex mapped_seqs_index(mapped_seqs_index_filepath, mapped_seqs_filepath);

  AllMappings all_mappings(mappings_filepath,
                           target_seqs_index,
                           MX_THRESHOLD_MIN,
                           MX_THRESHOLD_MAX,
                           mx_max_mapped_seqs_per_target_10kbp);

  serve(target_seqs_index,
        mapped_seqs_index,
        all_mappings,
        cbf_bytes,
        bf_bytes,
        BATCH_NAME_INPUT_PIPE,
        BATCH_TARGET_IDS_INPUT_READY_PIPE,
        TARGET_IDS_INPUT_PIPE,
        BFS_READY_PIPE,
        hash_num,
        k_values,
        subsample_max_mapped_seqs_per_target_10kbp);

  return 0;
}