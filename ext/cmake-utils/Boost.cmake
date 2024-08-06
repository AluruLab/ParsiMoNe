#### MPI
OPTION(USE_BOOST "Build with BOOST support" ON)
if (USE_BOOST)
    find_package(Boost REQUIRED 
        COMPONENTS system program_options filesystem log log_setup thread iostreams)
else(USE_BOOST)
    set(BOOST_FOUND 0)
endif(USE_BOOST)

if (Boost_FOUND)
    message(STATUS "Found Boost: ${Boost_VERSION_STRING}")
    message(STATUS "    headers: ${Boost_INCLUDE_DIRS}")
    message(STATUS "    lib dir: ${Boost_LIBRARY_DIRS}")
    message(STATUS "    libs:    ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIRS})
    set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})
else (Boost_FOUND)
    set(BOOST_DEFINE "")
endif (Boost_FOUND)
