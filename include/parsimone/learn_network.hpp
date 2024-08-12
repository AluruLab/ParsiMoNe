#ifndef LEARN_NETWORK_HPP
#define LEARN_NETWORK_HPP

#include <memory>
#include <mxx/comm.hpp>

#include "parsimone/ProgramOptions.hpp"
#include "common/DataReader.hpp"
#if __cplusplus >= 201703L // C++17 and later 
#include <string_view>
bool endsWith(std::string_view str, std::string_view suffix);
#else  // C++ 14 and earlier.
#include <string>
bool endsWith(const std::string& str, const char* suffix, unsigned suffixLen);
bool endsWith(const std::string& str, const char* suffix);
#endif

void learn_network(
  const ProgramOptions& options,
  const mxx::comm& comm,
  std::unique_ptr<DataReader<double>>&& reader
);

void learn_network(
  const ProgramOptions& options,
  const mxx::comm& comm,
  std::unique_ptr<DataReader<float>>&& reader
);

#endif // LEARN_NETWORK_HPP
