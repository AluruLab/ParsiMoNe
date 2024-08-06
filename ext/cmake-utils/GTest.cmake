
###### Google test library
# require "git submodule init" and "git submodule update"

# GTEST use of pthreads does not appear to mix well with openmp code, particularly clang.
# COMMENTED OUT 2/13/2023  set(gtest_compile_defs "-DGTEST_HAS_PTHREAD=0")
# COMMENTED OUT 2/13/2023  set(gtest_link_libs "gtest_main")

  
#------- include/load the Google Test framework
#  set_property(DIRECTORY ${PROJECT_SOURCE_DIR}/ext/gtest APPEND PROPERTY CMAKE_CXX_FLAGS -fPIC)
#  set_property(DIRECTORY ${PROJECT_SOURCE_DIR}/ext/gtest APPEND PROPERTY CMAKE_C_FLAGS -fPIC)
# COMMENTED OUT 2/13/2023  include_directories(${PROJECT_SOURCE_DIR}/ext/gtest/include)

# ADDED 2/13/2023 .  Version v.1.13.0.  
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

#------- enable the CMAKE testing framework
#  enable_testing()

# Set paths outside of the function so the function can
# be called from anywhere in the source tree
#   set path for output of testing binaries
#set(TEST_BINARY_OUTPUT_DIR ${CMAKE_BINARY_DIR}/test)
#   set path for output of test reports (XML format)
#set(TEST_XML_OUTPUT_DIR ${CMAKE_BINARY_DIR}/test_reports)

