include(CheckCXXSourceRuns)

set(CTZLL_CXX_CODE
    "int main() {
if( __builtin_ctzll(1 << 4) == 4)
    return 0;
else
    return 1;
}")

set(POPCOUNTLL_CXX_CODE
    "int main() {
if( __builtin_popcountll(0xFFFF) == 16)
    return 0;
else
    return 1;
}")


check_cxx_source_runs("${CTZLL_CXX_CODE}" BUILTIN_CTZLL_FOUND)
check_cxx_source_runs("${POPCOUNTLL_CXX_CODE}" BUILTIN_POPCOUNTLL_FOUND)


