#
# Armadillo
find_package(Armadillo REQUIRED)
if(Armadillo_FOUND)
    message(STATUS "Armadillo version:  ${ARMADILLO_VERSION_STRING} (${ARMADILLO_VERSION_NAME})")
    add_compile_definitions(ARMA_DONT_USE_WRAPPER)
    include_directories(${ARMADILLO_INCLUDE_DIR}) 
    set(EXTRA_LIBS ${EXTRA_LIBS} ${ARMADILLO_LIBRARIES})
endif()
#
# BLAS, LAPACK and LAPACKE
# BLAS
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    add_compile_definitions(ARMA_USE_BLAS)
    #target_link_libraries(${PROJECT_PYMOD} PUBLIC BLAS::BLAS)
    set(EXTRA_LIBS ${EXTRA_LIBS} BLAS::BLAS) 
endif()
# LAPACK
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    add_compile_definitions(ARMA_USE_LAPACK)
    # target_link_libraries(${PROJECT_PYMOD} PUBLIC LAPACK::LAPACK)
    set(EXTRA_LIBS ${EXTRA_LIBS} LAPACK::LAPACK) 
endif()

