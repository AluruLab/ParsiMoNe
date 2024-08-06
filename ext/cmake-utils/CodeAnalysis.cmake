###### CODE ANALYSIS FLAGS.  TODO: sanitizer flags for icc?  profiling flags for icc?
# enable_code_analysis - mutually exclusive options
#	enable_profiling - mutually exclusive options
#		enable_vtune_profiling
#   enable_google_profiling
#	enable_coverage
#	enable_sanitizer

OPTION(ENABLE_PROFILING "Enable Profiling" OFF)
OPTION(ENABLE_VTUNE_PROFILING "Turn on vtune specific profiling flags" OFF)
OPTION(ENABLE_COVERAGE "Enable Coverage" OFF)
#OPTION(ENABLE_GOOGLE_PROFILING "Turn on google perftools specific profiling flags" OFF)

# code analysis configuration
              
#### COMPILER based Profiling?
if(ENABLE_PROFILING)
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pg -g")
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pg -g")
else(ENABLE_PROFILING)
  # Remove unreferenced functions: function level linking
  if(NOT APPLE)
    add_definitions(-ffunction-sections)
  endif()
endif(ENABLE_PROFILING)

# vtunes profiling.
if(ENABLE_PROFILING)
  if (ENABLE_VTUNE_PROFILING)
    message(STATUS "VTUNE_ANALYSIS set to ${ENABLE_VTUNE_PROFILING}")
    ADD_LIBRARY(ittnotify STATIC IMPORTED)
    SET_TARGET_PROPERTIES(ittnotify PROPERTIES
        IMPORTED_LOCATION ${VTUNE_LIB})
  endif(ENABLE_VTUNE_PROFILING)

  # vtunes profiling.
  #message(STATUS "GPERFTOOLS_ANALYSIS set to ${ENABLE_GOOGLE_PROFILING}")
  #if (ENABLE_GOOGLE_PROFILING)
  #  ADD_LIBRARY(profiler SHARED IMPORTED)
  #	SET_TARGET_PROPERTIES(profiler PROPERTIES
  #	    IMPORTED_LOCATION /usr/lib/x86_64-linux-gnu)
  #endif(ENABLE_GOOGLE_PROFILING)
endif(ENABLE_PROFILING)


# coverage analysis
if (ENABLE_COVERAGE)
  # set flags for coverage test
  message(STATUS "Code coverage reporting enabled")
#      add_definitions(-fprofile-arcs -ftest-coverage -g)
  add_definitions(--coverage -g)
#      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -g")
#      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -g")
#      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -g")
#      set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -g")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage -g")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} --coverage -g")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage -g")
  set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} --coverage -g")
endif(ENABLE_COVERAGE)

