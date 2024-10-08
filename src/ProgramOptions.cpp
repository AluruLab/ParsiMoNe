/**
 * @file ProgramOptions.cpp
 * @brief Implementation of functionality for parsing command line options.
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
#include "parsimone/ProgramOptions.hpp"

#include <boost/filesystem.hpp>

#include <iostream>


namespace fs = boost::filesystem;

ProgramOptions::ProgramOptions(
) : m_desc("Parallel Construction of Module Networks"),
    m_logLevel(),
    m_logFile(),
    m_dataFile(),
    m_algoName(),
    m_outputDir(),
    m_configFile(),
    m_h5Path(),
    m_h5MatrixDataPath(),
    m_h5VarsDataPath(),
    m_h5ObsDataPath(),
    m_numVars(),
    m_numObs(),
    m_separator(),
    m_parallelRead(),
    m_colObs(),
    m_varNames(),
    m_obsIndices(),
    m_learnNetwork(),
    m_forceParallel(),
    m_hostNames(),
    m_warmupMPI()
{
  po::options_description basic("Basic options");
  basic.add_options()
    ("help,h", "Print this message")
    ("nvars,n", po::value<uint32_t>(&m_numVars), "Number of variables in the dataset")
    ("nobs,m", po::value<uint32_t>(&m_numObs), "Number of observations in the dataset")
    ("file,f", po::value<std::string>(&m_dataFile), "Name of the file from which dataset is to be read")
    ("readpar,r", po::bool_switch(&m_parallelRead)->default_value(false), "Read from the file in parallel")
    ("colobs,c", po::bool_switch(&m_colObs)->default_value(false), "The file contains observations in columns")
    ("separator,s", po::value<char>(&m_separator)->default_value(','), "Delimiting character in the file")
    ("varnames,v", po::bool_switch(&m_varNames)->default_value(false), "The file contains variable names")
    ("indices,i", po::bool_switch(&m_obsIndices)->default_value(false), "The file contains observation indices")
    ("algorithm,a", po::value<std::string>(&m_algoName)->default_value("lemontree"), "Name of the algorithm to be used")
    ("outdir,o", po::value<std::string>(&m_outputDir)->default_value("."), "Name of the directory to which the output files should be written")
    ("h5root", po::value<std::string>(&m_h5Path)->default_value("/"), "HDF5 Root Path for all data")
    ("h5matrix", po::value<std::string>(&m_h5MatrixDataPath)->default_value("matrix"), "HDF5 path to matrix data")
    ("h5obs", po::value<std::string>(&m_h5ObsDataPath)->default_value("col_attrs/CellID"), "HDF5 path to observations names")
    ("h5var", po::value<std::string>(&m_h5VarsDataPath)->default_value("row_attrs/Gene"), "HDF5 path to variable names")
    ;

  po::options_description advanced("Advanced options");
  advanced.add_options()
    ("config,g", po::value<std::string>(&m_configFile)->default_value(""), "JSON file with algorithm specific configurations")
    ("warmup,w", po::bool_switch(&m_warmupMPI)->default_value(false), "Warmup the MPI_Alltoall(v) functions before starting execution")
    ;

  po::options_description developer("Developer options");
  developer.add_options()
    ("parallel", po::bool_switch(&m_forceParallel)->default_value(false), "Use the parallel implementation even for p=1")
    ("hostnames", po::bool_switch(&m_hostNames)->default_value(false), "Print out the hostname for every process")
#ifdef LOGGING
    ("loglevel", po::value<std::string>(&m_logLevel)->default_value("error"), "Level of logging")
    ("logfile", po::value<std::string>(&m_logFile)->default_value(""), "File to which logs should be written")
#endif
    ;

  m_desc.add(basic).add(advanced).add(developer);
}

void
ProgramOptions::parse(
  int argc,
  char** argv
)
{
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, m_desc), vm);
  po::notify(vm);

  if ((argc == 1) || (vm.count("help") > 0)) {
    std::stringstream ss;
    ss << m_desc;
    throw po::error(ss.str());
  }
  if (!fs::exists(fs::path(m_dataFile))) {
    throw po::error("Couldn't find the data file");
  }
  if ((vm.count("nvars") == 0) || (vm.count("nobs") == 0)) {
    throw po::error("Dimensions of the data file should be provided using -n and -m");
  }
  if (m_configFile.empty()) {
    m_configFile = m_algoName + "_configs.json";
    std::cerr << "Using the default configuration file for the algorithm: " << m_configFile << std::endl;
  }
  if (!fs::exists(fs::path(m_configFile))) {
    throw po::error("Couldn't find the algorithm configuration file");
  }
}

uint32_t
ProgramOptions::numVars(
) const
{
  return m_numVars;
}

uint32_t
ProgramOptions::numObs(
) const
{
  return m_numObs;
}

const std::string&
ProgramOptions::dataFile(
) const
{
  return m_dataFile;
}

bool
ProgramOptions::parallelRead(
) const
{
  return m_parallelRead;
}

bool
ProgramOptions::colObs(
) const
{
  return m_colObs;
}

bool
ProgramOptions::varNames(
) const
{
  return m_varNames;
}

bool
ProgramOptions::obsIndices(
) const
{
  return m_obsIndices;
}

char
ProgramOptions::separator(
) const
{
  return m_separator;
}

const std::string&
ProgramOptions::algoName(
) const
{
  return m_algoName;
}

bool
ProgramOptions::learnNetwork(
) const
{
  return m_learnNetwork;
}

const std::string&
ProgramOptions::outputDir(
) const
{
  return m_outputDir;
}

const std::string&
ProgramOptions::configFile(
) const
{
  return m_configFile;
}

bool
ProgramOptions::forceParallel(
) const
{
  return m_forceParallel;
}

bool
ProgramOptions::hostNames(
) const
{
  return m_hostNames;
}

bool
ProgramOptions::warmupMPI(
) const
{
  return m_warmupMPI;
}

const std::string&
ProgramOptions::logLevel(
) const
{
  return m_logLevel;
}

const std::string&
ProgramOptions::logFile(
) const
{
  return m_logFile;
}

const std::string&
ProgramOptions::h5root(
) const
{
  return m_h5Path;
}

const std::string&
ProgramOptions::h5matrixPath(
) const
{
  return m_h5MatrixDataPath;
}

const std::string&
ProgramOptions::h5obsPath(
) const
{
  return m_h5VarsDataPath;
}

const std::string&
ProgramOptions::h5varPath(
) const
{
  return m_h5VarsDataPath;
}

ProgramOptions::~ProgramOptions(
)
{
}
