#include "utils.hpp"

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/data_stream.hpp"
#include "btllib/nthash.hpp"
#include "btllib/status.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static const pid_t INIT_PID = 1;
static const unsigned PARENT_QUERY_PERIOD = 1; // seconds
static const auto FIFO_FLAGS = S_IRUSR | S_IWUSR;

std::vector<size_t>
get_random_indices(const size_t total_size, const size_t count)
{
  btllib::check_error(count > total_size,
                      FN_NAME + ": count cannot be larger than total_size.");

  static std::random_device dev;
  static std::mt19937 rng(dev());

  std::vector<size_t> all_indices(total_size);
  std::iota(all_indices.begin(), all_indices.end(), 0);
  std::shuffle(all_indices.begin(), all_indices.end(), rng);

  using advance_type = typename std::iterator_traits<
    decltype(all_indices.begin())>::difference_type;

  decltype(all_indices) indices(
    all_indices.begin(), std::next(all_indices.begin(), advance_type(count)));

  return indices;
}

void
wait_till_parent_ends()
{
  while (getppid() != INIT_PID) {
    std::this_thread::sleep_for(std::chrono::seconds(PARENT_QUERY_PERIOD));
  }
}

void
bind_to_parent()
{
  (new std::thread([]() {
    wait_till_parent_ends();
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }))
    ->detach();
}

void
make_pipe(const std::string& pipepath)
{
  const auto ret = mkfifo(pipepath.c_str(), FIFO_FLAGS);
  btllib::check_error(ret != 0,
                      FN_NAME + ": mkfifo failed: " + btllib::get_strerror());
}

// Pipe messages are newline-terminated lines, and readers stop at the last line
// they expect instead of reading to end of file. On macOS, the final
// end-of-file read() on a named pipe can block forever after the writer has
// closed it, so a reader must never depend on that read. Raw read() is used
// because stdio and iostream buffering keep reading until the buffer is full or
// end of file is reached.
std::vector<std::string>
read_pipe_lines(const std::string& pipepath, const std::string& end_line)
{
  const int fd = open(pipepath.c_str(), O_RDONLY | O_CLOEXEC);
  btllib::check_error(fd < 0,
                      FN_NAME + ": failed to open " + pipepath + ": " +
                        btllib::get_strerror());
  std::vector<std::string> lines;
  std::string pending;
  std::vector<char> buf(4096); // NOLINT(readability-magic-numbers,cppcoreguidelines-avoid-magic-numbers)
  bool done = false;
  while (!done) {
    const auto n = read(fd, buf.data(), buf.size());
    if (n < 0 && errno == EINTR) {
      continue;
    }
    btllib::check_error(n < 0,
                        FN_NAME + ": failed to read " + pipepath + ": " +
                          btllib::get_strerror());
    if (n == 0) { // Writer closed without sending end_line.
      break;
    }
    pending.append(buf.data(), n);
    size_t start = 0;
    size_t newline = 0;
    while ((newline = pending.find('\n', start)) != std::string::npos) {
      auto line = pending.substr(start, newline - start);
      start = newline + 1;
      if (end_line.empty()) {
        lines.push_back(std::move(line));
        done = true;
        break;
      }
      if (line == end_line) {
        done = true;
        break;
      }
      lines.push_back(std::move(line));
    }
    pending.erase(0, start);
  }
  if (!done && !pending.empty()) {
    lines.push_back(pending);
  }
  close(fd);
  return lines;
}

std::string
read_pipe(const std::string& pipepath)
{
  const auto lines = read_pipe_lines(pipepath, "");
  if (lines.empty()) {
    return "";
  }
  // Same token as the former `input >> str`: the first whitespace-free word.
  std::istringstream ss(lines.front());
  std::string str;
  ss >> str;
  return str;
}

void
confirm_pipe(const std::string& pipepath)
{
  std::ofstream confirm(pipepath);
  confirm << "1" << std::endl;
}

void
fill_bfs(const char* seq,
         const size_t seq_len,
         const unsigned hash_num,
         const std::vector<unsigned>& k_values,
         const unsigned kmer_threshold,
         std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
         std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs)
{
  btllib::check_error(kmer_threshold < 4,
                      FN_NAME + ": kmer_threshold must be "
                                "greater than or equal to 4.");
  unsigned adjusted_kmer_threshold = kmer_threshold - 2;
  for (size_t i = 0; i < k_values.size(); i++) {
    const auto k = k_values[i];
    const auto& cbf = cbfs[i];
    const auto& bf = bfs[i];
    btllib::NtHash nthash(seq, seq_len, hash_num, k);
    while (nthash.roll()) {
      if (cbf->insert_thresh_contains(nthash.hashes(),
                                      adjusted_kmer_threshold) >=
          adjusted_kmer_threshold) {
        bf->insert(nthash.hashes());
      }
    }
    adjusted_kmer_threshold++;
  }
}
