#### native hardware architecture
OPTION(USE_SIMD "Enable SIMD instructions, if available on hardware. (-march=native)" ON)
if (USE_SIMD)
    # add_compile_options(-march=native)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp-simd -ftree-vectorize")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd -ftree-vectorize")    
        # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fabi-version=0")
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fabi-version=0")    
    endif()
    
    add_definitions(-DUSE_SIMD)
endif(USE_SIMD)

# TODO autovectorization
