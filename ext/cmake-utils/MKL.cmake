OPTION(USE_MKL "whether lightpcc should use MKL" ON)
if (USE_MKL)
    # single dynamic library option, based on https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
    set(MKL_INCLUDE_DIRS "$ENV{MKLROOT}/include")
    set(MKL_LIB_DIRS "$ENV{MKLROOT}/lib/intel64")
    #set(MKL_LIBS "mkl_rt")                     # implicit version.
    #set(MKL_COMPILE_FLAGS "-mkl=parallel")     # implicit version.
    set(MKL_COMPILE_FLAGS "")                   # explicit version
    set(MKL_LINK_FLAGS "")
    set(MKL_DEFINES "-DMKL_ILP64")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(MKL_LIBS "mkl_intel_ilp64;mkl_gnu_thread;mkl_core")          # explicit version
        set(MKL_COMPILE_FLAGS "-m64")
        set(MKL_LINK_FLAGS "-Wl,--no-as-needed")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(MKL_LIBS "mkl_intel_ilp64;mkl_intel_thread;mkl_core")       # explicit version
    endif()
endif (USE_MKL)
