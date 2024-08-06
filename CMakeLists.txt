cmake_minimum_required(VERSION 3.14)
include(CheckIncludeFile)
include(CheckSymbolExists)

# project settings
project(parsimone C CXX)

set(PARSIMONE_VERSION_MAJOR 1)
set(PARSIMONE_VERSION_MINOR 0)
set(PARSIMONE_VERSION_PATCH 0)
set(PARSIMONE_VERSION_TWEAK 0)
set(PARSIMONE_APP parsimone)
set(PARSIMONE_TEST_APP parsimone_test)

set(CMAKE_VERBOSE_MAKEFILE ON)

# Requires at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Initialize CXXFLAGS and CFLAGS.
set(BASE_OPT_FLAGS "-ffast-math -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${BASE_OPT_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${BASE_OPT_FLAGS}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fPIC")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${BASE_OPT_FLAGS}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${BASE_OPT_FLAGS}")

# Specify build type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
  SET(ENABLE_LOGGING ON)
ENDIF()

# Building Options 
OPTION(BUILD_TESTS "whether Tests should be built" OFF)
OPTION(ENABLE_LOGGING "whether Logging should be built" OFF)
OPTION(ENABLE_PROFILING "whether Profiling should be built" OFF)
OPTION(ENABLE_SANITIZER "whether Sanitizer should be enabled" OFF)
OPTION(ENABLE_TIMING "whether Timing should be enabled" OFF)
OPTION(PRINT_VEC_REPORT "whether Vectorization report should be enabled" OFF)

# Warning Options
set(COMPILER_WARNINGS ON)
set(COMPILER_EXTRA_WARNINGS ON)
set(COMPILER_WARNINGS_ENABLE_EXTRA ON)
set(COMPILER_WARNINGS_ENABLE_CONVERSIONS OFF)
set(COMPILER_WARNINGS_ENABLE_SUGGESTIONS OFF)

# Sanitizer or Profiling
if (ENABLE_SANITIZER)
    set(SANITIZER_STYLE address leak undefined)
endif()

# External directories
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(CMAKE_MODULE_DIR ${PROJECT_SOURCE_DIR}/ext/cmake-utils)

# Save libraries and executables in the same place
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE PATH "Output directory for applications" )

# Add the binary tree to the search path for include files
include_directories("${PROJECT_SOURCE_DIR}/include")
# Add external dependencis to the search path
include_directories("${EXT_PROJECTS_DIR}")
include_directories("${EXT_PROJECTS_DIR}/mxx/include")
include_directories("${EXT_PROJECTS_DIR}/mxx/ext/")

# mm_malloc.h header should exist and _mm_malloc function should be available
check_include_file(mm_malloc.h MM_MALLOC_HDR_EXISTS)
check_symbol_exists(_mm_malloc "mm_malloc.h" MM_MALLOC_FN_EXISTS)
if(NOT (MM_MALLOC_HDR_EXISTS AND MM_MALLOC_FN_EXISTS))
    message(FATAL_ERROR "_mm_malloc is required by the bitvector counter from SABNAtk")
endif()

include("${CMAKE_MODULE_DIR}/ExtraWarnings.cmake")
include("${CMAKE_MODULE_DIR}/CompilerProfiling.cmake")
# include("${CMAKE_MODULE_DIR}/OpenMP.cmake")
include("${CMAKE_MODULE_DIR}/MPI.cmake")
include("${CMAKE_MODULE_DIR}/SIMD.cmake")
include("${CMAKE_MODULE_DIR}/Builtins.cmake")
include("${CMAKE_MODULE_DIR}/AVX.cmake")
include("${CMAKE_MODULE_DIR}/SSE.cmake")
include("${CMAKE_MODULE_DIR}/Boost.cmake")
include("${CMAKE_MODULE_DIR}/Arma.cmake")
# include("${CMAKE_MODULE_DIR}/HDF5.cmake")
include("${CMAKE_MODULE_DIR}/CodeAnalysis.cmake")
include("${CMAKE_MODULE_DIR}/Sanitizer.cmake")
#
# trng4
set(TRNG_ENABLE_TESTS OFF)
set(TRNG_ENABLE_EXAMPLES OFF)
add_subdirectory(ext/trng4)
set(EXTRA_LIBS ${EXTRA_LIBS} trng4)

# Check for builitings
if(NOT (BUILTIN_CTZLL_FOUND AND BUILTIN_POPCOUNTLL_FOUND))
    message(FATAL_ERROR "__builtin_ctzll and __builtin_popcountll required by bit_util.hpp")
endif()

# Definitions, Flags and Libraries
set(app_compile_defs "-DVERBOSE")
set(app_compile_flags "")
set(app_link_flags "")
set(app_link_libs Boost::system Boost::program_options Boost::filesystem Boost::iostreams)

# One of AVX2, AVX512 BW, SSE 4.1 or SSE 4.2 should be supported
if(NOT (AVX2_SUPPORTED OR AVX512BW_SUPPORTED OR SSE41_SUPPORTED OR SSE42_SUPPORTED))
    message(WARNING "No Vectorization Support: Performance of the SABNAtk counters will be impacted")
endif()

if(AVX2_SUPPORTED)
    set(app_compile_defs "${app_compile_defs};-DHAVE_AVX2_INSTRUCTIONS")
endif()

if(AVX512BW_SUPPORTED)
    set(app_compile_defs "${app_compile_defs};-DHAVE_AVX512BW_INSTRUCTIONS")
endif()

# Vectorization report
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(app_compile_defs "${app_compile_defs};-DHAVE_INTEL_ICC")
    if (PRINT_VEC_REPORT)
        set(app_compile_flags "${app_compile_flags};-qopt-report=5;-qopt-report-phase=vec;-qopt-report-file=stdout")
    endif(PRINT_VEC_REPORT)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (PRINT_VEC_REPORT)
        set(app_compile_flags "${app_compile_flags};-ftree-vectorizer-verbose=4")
    endif(PRINT_VEC_REPORT)
endif()

# If Logging enabled, add dependencies to compile definitions and link flags
if(ENABLE_LOGGING)
    set(app_compile_defs "${app_compile_defs};-DLOGGING")
    set(app_link_libs ${app_link_libs} Boost::log Boost::log_setup 
        Boost::thread)
endif(ENABLE_LOGGING)

if(ENABLE_TIMING)
  set(app_compile_defs "${app_compile_defs};-DTIMER")
endif(ENABLE_TIMING)

# If sanitizer enabled, add sanitizer compile flags and link flags
if(ENABLE_SANITIZER)
    set(app_compile_flags "${app_compile_flags};${SANITIZER_COMPILE_FLAGS}")
    set(app_link_flags "${app_link_flags};${SANITIZER_LINK_FLAGS}")
endif()

add_executable(${PARSIMONE_APP} src/parsimone.cpp src/ProgramOptions.cpp)
foreach (cdef IN LISTS app_compile_defs)
    target_compile_definitions(${PARSIMONE_APP} PRIVATE ${cdef})
endforeach(cdef)
foreach (cflgs IN LISTS app_compile_flags)
    target_compile_options(${PARSIMONE_APP} PRIVATE ${cflgs})
endforeach(cflgs)
target_link_libraries(${PARSIMONE_APP} ${app_link_flags} ${app_link_libs} ${EXTRA_LIBS})

# Building tests
if (BUILD_TESTS)
    include_directories("${PROJECT_SOURCE_DIR}/test")
    include_directories("${EXT_PROJECTS_DIR}/mxx/gtest")
    cmake_policy(SET CMP0135 NEW)
    include("${CMAKE_MODULE_DIR}/GTest.cmake")
    # add_subdirectory("ext/googletest")
    add_executable(${PARSIMONE_TEST_APP} test/test.cpp)
    foreach (cdef IN LISTS app_compile_defs)
        target_compile_definitions(${PARSIMONE_TEST_APP} PRIVATE ${cdef})
    endforeach(cdef)
    foreach (cflgs IN LISTS app_compile_flags)
        target_compile_options(${PARSIMONE_TEST_APP} PRIVATE ${cflgs})
    endforeach(cflgs)
    target_link_libraries(${PARSIMONE_TEST_APP} ${app_link_flags} ${app_link_libs} ${EXTRA_LIBS} GTest::gtest_main)
endif()
