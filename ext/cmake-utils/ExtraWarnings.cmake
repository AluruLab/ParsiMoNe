###### Extra compiler flags and prettifying.
include(CMakeDependentOption)
#
OPTION(COMPILER_EXTRA_WARNINGS "EXTRA COMPILATION OPTIONS" OFF)
CMAKE_DEPENDENT_OPTION(COMPILER_WARNINGS_ENABLE_EXTRA "Enable compiler to generate extra warnings?" OFF "COMPILER_EXTRA_WARNINGS" OFF)
CMAKE_DEPENDENT_OPTION(COMPILER_WARNINGS_ENABLE_CONVERSIONS "Enable compiler to generate extra warnings for type conversions?" OFF "COMPILER_WARNINGS" OFF)
CMAKE_DEPENDENT_OPTION(COMPILER_WARNINGS_ENABLE_SUGGESTIONS "Enable compiler to generate suggestions for polymorphic type analysis?" OFF "COMPILER_WARNINGS" OFF)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  SET(EXTRA_WARNING_FLAGS "-Wextra -Wno-unused-parameter -Wcheck")
  SET(TYPE_CONVERSION_WARNING_FLAGS "-Wconversion")
  SET(SUGGESTION_WARNING_FLAGS "-Wremarks")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # -Wzero-as-null-pointer-constant : causes a lot of errors in system headers.
  # -Wfloat-equal : most of comparisons are not with float.
  SET(EXTRA_WARNING_FLAGS "-Wextra -Wno-unused-parameter") #-Wundef 
  SET(TYPE_CONVERSION_WARNING_FLAGS "-Wdouble-promotion -Wconversion -Wsign-conversion -Wcast-qual -Wuseless-cast")
  SET(SUGGESTION_WARNING_FLAGS "-Wsuggest-override -Wsuggest-final-types -Wsuggest-final-methods")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET(EXTRA_WARNING_FLAGS "-Wextra -Wno-unused-parameter" CACHE INTERNAL "extra compiler warning flags")
  SET(TYPE_CONVERSION_WARNING_FLAGS "-W" CACHE INTERNAL "compiler flags to check type conversions")
  SET(SUGGESTION_WARNING_FLAGS CACHE INTERNAL "compiler flags that suggest keywords for better type resolutions")
endif()

if(COMPILER_WARNINGS_ENABLE_EXTRA)
    add_definitions(${EXTRA_WARNING_FLAGS})
endif(COMPILER_WARNINGS_ENABLE_EXTRA)

if(COMPILER_WARNINGS_ENABLE_CONVERSIONS)
    add_definitions(${TYPE_CONVERSION_WARNING_FLAGS})
endif(COMPILER_WARNINGS_ENABLE_CONVERSIONS)

if(COMPILER_WARNINGS_ENABLE_SUGGESTIONS)
    add_definitions(${SUGGESTION_WARNING_FLAGS})
endif(COMPILER_WARNINGS_ENABLE_SUGGESTIONS)

