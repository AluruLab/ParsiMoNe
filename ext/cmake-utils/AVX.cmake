# Include the required headers
include(CheckCXXSourceCompiles)

# From : https://trycatchdebug.net/news/1127634/avx2-avx-512-detection-in-cmake
# Define the source code snippet for AVX2 detection
set(AVX2_CXX_CODE
    "#include <immintrin.h>
int main() {
#if __AVX2__
return 0;
#else
#error \"AVX2 is not supported\"
#endif
}")

# Define the source code snippet for AVX-512F detection
set(AVX512F_CXX_CODE
    "#include <immintrin.h>
int main() {
#if __AVX512F__
return 0;
#else
#error \"AVX-512 F is not supported\"
#endif
}")

# Define the source code snippet for AVX-512F detection
set(AVX512BW_CXX_CODE
    "#include <immintrin.h>
int main() {
#if __AVX512BW__
return 0;
#else
#error \"AVX-512 BW is not supported\"
#endif
}")

# Check for AVX2 support
check_cxx_source_compiles("${AVX2_CXX_CODE}" AVX2_SUPPORTED)
# Check for AVX-512F support
check_cxx_source_compiles("${AVX512F_CXX_CODE}" AVX512F_SUPPORTED)
# Check for AVX-512BW support
check_cxx_source_compiles("${AVX512BW_CXX_CODE}" AVX512BW_SUPPORTED)

#
if(AVX2_SUPPORTED)
  add_definitions(-mavx2)
endif()

if(AVX512F_SUPPORTED)
  add_definitions(-mavx512f)
endif()

if(AVX512BW_SUPPORTED)
  add_definitions(-mavx512bw)
endif()
