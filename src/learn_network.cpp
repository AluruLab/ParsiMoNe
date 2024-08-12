/**
 * @file parsimone.cpp
 * @brief The implementation of the main function for ParsiMoNe,
 *        and other functions that drive the program execution.
 * @author Ankit Srivastava <asrivast@gatech.edu>
 *
 * Copyright 2020 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "parsimone/RawData.hpp"
#include "parsimone/Genomica.hpp"
#include "parsimone/LemonTree.hpp"
#include "parsimone/ProgramOptions.hpp"

#include "common/UintSet.hpp"
#include "utils/Timer.hpp"

#include <boost/asio/ip/host_name.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>
#include <memory>


#if __cplusplus >= 201703L // C++17 and later 
#include <string_view>
bool endsWith(std::string_view str, std::string_view suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}
#else  // C++ 14 and earlier.
bool endsWith(const std::string& str, const char* suffix, unsigned suffixLen)
{
    return str.size() >= suffixLen && 0 == str.compare(str.size()-suffixLen, suffixLen, suffix, suffixLen);
}

bool endsWith(const std::string& str, const char* suffix)
{
    return endsWith(str, suffix, std::string::traits_type::length(suffix));
}
#endif


/**
 * @brief Gets a pointer to the object of the required module network learning algorithm.
 *
 * @tparam Var Type of the variables (expected to be an integral type).
 * @tparam Set Type of set container.
 * @tparam Data Type of the object which is used for querying data.
 * @param algoName The name of the algorithm.
 * @param data The object which is used for querying data.
 *
 * @return unique_ptr to the object of the given algorithm.
 *         The unique_ptr points to a nullptr if the algorithm is not found.
 */
template <typename Var, typename Set, typename Data>
std::unique_ptr<ModuleNetworkLearning<Data, Var, Set>>
getAlgorithm(
  const std::string& algoName,
  const mxx::comm& comm,
  const Data& data
)
{
  std::stringstream ss;
  if (algoName.compare("lemontree") == 0) {
    return std::make_unique<LemonTree<Data, Var, Set>>(comm, data);
  }
  ss << "lemontree";
  if (algoName.compare("genomica") == 0) {
    return std::make_unique<Genomica<Data, Var, Set>>(comm, data);
  }
  ss << ",genomica";
  throw std::runtime_error("Requested algorithm not found. Supported algorithms are: {" + ss.str() + "}");
  return std::unique_ptr<ModuleNetworkLearning<Data, Var, Set>>();
}

pt::ptree
readConfigs(
  const std::string& configFile,
  const mxx::comm& comm
)
{
  TIMER_DECLARE(tConfigs);
  std::string configStr;
  std::stringstream buffer;
  if (comm.is_first()) {
    std::ifstream cf(configFile);
    buffer << cf.rdbuf();
    configStr = buffer.str();
  }
  mxx::bcast(configStr, 0, comm);
  if (!comm.is_first()) {
    buffer << configStr;
  }
  namespace pt = boost::property_tree;
  pt::ptree configs;
  pt::read_json(buffer, configs);
  if (comm.is_first()) {
    TIMER_ELAPSED("Time taken in reading the configs: ", tConfigs);
  }
  return configs;
}

/**
 * @brief Learns the module network with the given parameters
 *        and writes it to the given file.
 *
 * @tparam Var Type of the variables (expected to be an integral type).
 * @param options Program options provider.
 */
template <typename Var, typename Size, typename Data>
void
learnNetwork(
  const ProgramOptions& options,
  const mxx::comm& comm,
  const Data& data
)
{
  auto algo = getAlgorithm<Var, UintSet<Var, Size>>(options.algoName(), comm, data);
  auto configs = readConfigs(options.configFile(), comm);
  if (comm.is_first()) {
    namespace fs = boost::filesystem;
    if (!fs::is_directory(options.outputDir())) {
      if (!fs::create_directories(options.outputDir())) {
        throw po::error("Output directory doesn't exist and could not be created");
      }
    }
    fs::copy_file(fs::path(options.configFile()), options.outputDir() + "/configs.json", fs::copy_options::overwrite_existing);
  }
  comm.barrier();
  TIMER_DECLARE(tNetwork);
  algo->learnNetwork((comm.size() > 1) || options.forceParallel(), configs, options.outputDir());
  comm.barrier();
  if (comm.is_first()) {
    TIMER_ELAPSED("Time taken in getting the network: ", tNetwork);
  }
}

/**
 * @brief Learns the module network with the given parameters
 *        and writes it to the given file.
 *
 * @tparam FileType Type of the file to be read.
 * @param options Program options provider.
 * @param reader File data reader.
 */
template <template <typename> class Reader, typename DataType>
void
learnNetwork(
  const ProgramOptions& options,
  const mxx::comm& comm,
  std::unique_ptr<Reader<DataType>>&& reader
)
{
  auto n = options.numVars();
  auto m = options.numObs();
  auto s = std::max(n, m);
  if ((s - 1) <= UintSet<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>::capacity()) {
    RawData<DataType, uint8_t> data(reader->data(), reader->varNames(), static_cast<uint8_t>(n), static_cast<uint8_t>(m));
    learnNetwork<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>::capacity()) {
    RawData<DataType, uint8_t> data(reader->data(), reader->varNames(), static_cast<uint8_t>(n), static_cast<uint8_t>(m));
    learnNetwork<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint8_t>::capacity()) {
    RawData<DataType, uint8_t> data(reader->data(), reader->varNames(), static_cast<uint8_t>(n), static_cast<uint8_t>(m));
    learnNetwork<uint8_t, std::integral_constant<int, maxSize<uint8_t>()>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, maxSize<uint16_t>()>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, maxSize<uint16_t>()>>(options, comm, data);
  }
  else {
    throw std::runtime_error("The given number of variables and observations is not supported.");
  }
}

#include "common/DataReader.hpp"

void learn_network(
  const ProgramOptions& options,
  const mxx::comm& comm,
  std::unique_ptr<DataReader<double>>&& reader
){
    learnNetwork(options, comm, std::move(reader));
}

void learn_network(
  const ProgramOptions& options,
  const mxx::comm& comm,
  std::unique_ptr<DataReader<float>>&& reader
){
    learnNetwork(options, comm, std::move(reader));
}
