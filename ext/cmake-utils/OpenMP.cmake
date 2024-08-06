#### OpenMP
#include(FindOpenMP)
find_package(OpenMP)

if (OPENMP_FOUND)
    message(STATUS "Found OpenMP")
else(OPENMP_FOUND)
	message(STATUS "NO OpenMP.  check compiler version now.")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    EXECUTE_PROCESS( COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE clang_full_version_string )
    string (REGEX REPLACE ".*clang version ([0-9]+\\.[0-9]+).*" "\\1" CLANG_VERSION ${clang_full_version_string})
  
   	if (CLANG_VERSION VERSION_GREATER 3.7 OR CLANG_VERSION VERSION_EQUAL 3.7)
	  	message(STATUS "Found OpenMP for CLANG ${CLANG_VERSION}")
	  	set(OPENMP_FOUND 1)
	  	set(OpenMP_C_FLAGS "-fopenmp=libomp")
	  	set(OpenMP_CXX_FLAGS "-fopenmp=libomp")
	  	
   	    set(CLANG_OPENMP_HOME "${CLANG_COMPILER_DIR}/.." CACHE PATH "Path to Clang OpenMP root directory.  Ideally, same as clang home.")
   	    include_directories(${CLANG_OPENMP_HOME}/include)  # clang uses gcc headers, but can't find all gcc headers by itself.  this is for omp.h
   	    link_directories(${CLANG_OPENMP_HOME}/lib)
   	
	  else()
  		set(OPNEMP_FOUND 0)
#	    message(FATAL_ERROR "${PROJECT_NAME} requires clang 3.7 or greater for OpenMP support.")
  	endif()
  endif() 
endif(OPENMP_FOUND)

CMAKE_DEPENDENT_OPTION(USE_OPENMP "Build with OpenMP support" ON
                        "OPENMP_FOUND" OFF)
#variable_watch(USE_OPENMP update_logging_engine)

if (USE_OPENMP)
    # add OpenMP flags to compiler flags
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        # FindOpenMP seems to want to set the flag to -fopenmp, which does not match intel documentation.
        set(EXTRA_LIBS "${EXTRA_LIBS};iomp5")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(EXTRA_LIBS "${EXTRA_LIBS};gomp")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
    add_definitions(-DUSE_OPENMP)
else (USE_OPENMP)
    set(OPENMP_DEFINE "")
endif (USE_OPENMP)


# if OMP_DEBUGGING is enabled, default(none) is set in the omp pragmas.   this means that only NO_LOG and PRINTF logging are compatible.
CMAKE_DEPENDENT_OPTION(OPENMP_STRICT_SCOPING "Enable OpenMP debugging (turns on default(none) in omp pragma)" ON
                        "OMP_DEBUGGING" OFF)
mark_as_advanced(OPENMP_STRICT_SCOPING)
if (OPENMP_STRICT_SCOPING)
    add_definitions(-DOMP_DEBUG)
endif(OPENMP_STRICT_SCOPING)
