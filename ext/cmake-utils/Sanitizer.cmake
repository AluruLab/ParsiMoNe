# Option for Sanitizer
OPTION(ENABLE_SANITIZER "Enable Santization" OFF)

if(NOT (DEFINED SANITIZER_STYLE))
    set(SANITIZER_STYLE address leak undefined)
endif()

if (ENABLE_SANITIZER)
  # set flags for coverage test
  message(STATUS "Sanitizer enabled with ${SANITIZER_STYLE}")
  # probably should be specific based on sanitizer style.  no-omit-frame-pointer is for address.
  set(SANITIZER_COMPILE_FLAGS "-fPIE;-fno-omit-frame-pointer;-g")
  foreach(X IN LISTS SANITIZER_STYLE)
    set(SANITIZER_COMPILE_FLAGS "-fsanitize=${X};${SANITIZER_COMPILE_FLAGS}")
  endforeach()

  set(SANITIZER_LINK_FLAGS "-fno-omit-frame-pointer;-pie")
  foreach(X IN LISTS SANITIZER_STYLE)
    set(SANITIZER_LINK_FLAGS "-fsanitize=${X};${SANITIZER_LINK_FLAGS}")
  endforeach()

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if ("thread" IN_LIST SANITIZER_STYLE)
      set(SANITIZER_EXTRA_FLAGS "-static-libtsan")
      set(SANITIZER_LINK_FLAGS "${SANITIZER_LINK_FLAGS};-static-libtsan")
    else()
        set(SANITIZER_EXTRA_FLAGS "")
    endif()
  else()
    set(SANITIZER_EXTRA_FLAGS "")	
  endif()
  

  # NOTE: if using thread sanitizer, please use g++ 4.9 and later, compiled with --disable-linux-futex (for libgomp)
  #  else lots of false positives from OMP will be reported by tsan
else()
  set(SANITIZER_COMPILE_FLAGS "")
  set(SANITIZER_LINK_FLAGS "")
endif()
