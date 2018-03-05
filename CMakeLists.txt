cmake_minimum_required(VERSION 3.0)
project(ucr-suite-cpp C CXX)

# silence policy CMP0042 warning by enabling RPATH explicitly
if(APPLE AND NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH true)
endif()

# Appends `str` to the variable named `var` with a whitespace as separator.
# Suppresses a leading whitespace if the variable is empty and does nothing if
# `str` is empty.
function(build_string var str)
  if(NOT str STREQUAL "")
    if("${${var}}" STREQUAL "")
      set("${var}" "${str}" PARENT_SCOPE)
    else()
      set("${var}" "${${var}} ${str}" PARENT_SCOPE)
    endif()
  endif()
endfunction(build_string)

# prohibit in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed. Please use "
                        "./configure to choose a build directory and "
                        "initialize the build configuration.")
endif ()

if(UCR_CXX_WARNINGS_AS_ERRORS)
  build_string("EXTRA_FLAGS" "-Werror")
endif()

if(UCR_MORE_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(WFLAGS "-Weverything -Wno-c++98-compat -Wno-padded "
               "-Wno-documentation-unknown-command -Wno-exit-time-destructors "
               "-Wno-global-constructors -Wno-missing-prototypes "
               "-Wno-c++98-compat-pedantic -Wno-unused-member-function "
               "-Wno-unused-const-variable -Wno-switch-enum "
               "-Wno-missing-noreturn -Wno-covered-switch-default")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(WFLAGS "-Waddress -Wall -Warray-bounds "
               "-Wattributes -Wbuiltin-macro-redefined -Wcast-align "
               "-Wcast-qual -Wchar-subscripts -Wclobbered -Wcomment "
               "-Wconversion -Wconversion-null -Wcoverage-mismatch "
               "-Wcpp -Wdelete-non-virtual-dtor -Wdeprecated "
               "-Wdeprecated-declarations -Wdiv-by-zero -Wdouble-promotion "
               "-Wempty-body -Wendif-labels -Wenum-compare -Wextra "
               "-Wfloat-equal -Wformat -Wfree-nonheap-object "
               "-Wignored-qualifiers -Winit-self "
               "-Winline -Wint-to-pointer-cast -Winvalid-memory-model "
               "-Winvalid-offsetof -Wlogical-op -Wmain -Wmaybe-uninitialized "
               "-Wmissing-braces -Wmultichar "
               "-Wnarrowing -Wnoexcept -Wnon-template-friend "
               "-Wnon-virtual-dtor -Wnonnull -Woverflow "
               "-Woverlength-strings -Wparentheses "
               "-Wpmf-conversions -Wpointer-arith -Wreorder "
               "-Wreturn-type -Wsequence-point -Wshadow "
               "-Wsign-compare -Wswitch -Wtype-limits -Wundef "
               "-Wuninitialized -Wunused -Wvla -Wwrite-strings")
  endif()
  # convert CMake list to a single string, erasing the ";" separators
  string(REPLACE ";" "" WFLAGS_STR ${WFLAGS})
  build_string("EXTRA_FLAGS" "${WFLAGS_STR}")
endif()

if(NOT UCR_NO_AUTO_LIBCPP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CXXFLAGS_BACKUP "${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "-std=c++1z -stdlib=libc++")
  try_run(ProgramResult
          CompilationSucceeded
          "${CMAKE_CURRENT_BINARY_DIR}"
          "${CMAKE_CURRENT_SOURCE_DIR}/cmake/get_compiler_version.cpp"
          RUN_OUTPUT_VARIABLE CompilerVersion)
  if(NOT CompilationSucceeded OR NOT ProgramResult EQUAL 0)
    message(STATUS "Use clang with GCC' libstdc++")
  else()
    message(STATUS "Automatically added '-stdlib=libc++' flag "
      "(UCR_NO_AUTO_LIBCPP not defined)")
  build_string("EXTRA_FLAGS" "-stdlib=libc++")
  endif()
  # restore CXX flags
  set(CMAKE_CXX_FLAGS "${CXXFLAGS_BACKUP}")
endif()

if(UCR_ENABLE_ADDRESS_SANITIZER)
  # check whether address sanitizer is available
  set(CXXFLAGS_BACKUP "${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "-fsanitize=address -fno-omit-frame-pointer")
  try_run(ProgramResult
          CompilationSucceeded
          "${CMAKE_CURRENT_BINARY_DIR}"
          "${CMAKE_CURRENT_SOURCE_DIR}/cmake/get_compiler_version.cpp")
  if(NOT CompilationSucceeded)
    message(WARNING "Address Sanitizer is not available on selected compiler")
  else()
    message(STATUS "Enable Address Sanitizer")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()
  # restore CXX flags
  set(CMAKE_CXX_FLAGS "${CXXFLAGS_BACKUP}")
endif(UCR_ENABLE_ADDRESS_SANITIZER)

# check if the user provided CXXFLAGS, set defaults otherwise
if(NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS                   "-std=c++1z -Wextra -Wall -pedantic ${EXTRA_FLAGS}")
endif()
if(NOT CMAKE_CXX_FLAGS_DEBUG)
  set(CMAKE_CXX_FLAGS_DEBUG             "-O0 -g")
endif()
if(NOT CMAKE_CXX_FLAGS_MINSIZEREL)
  set(CMAKE_CXX_FLAGS_MINSIZEREL        "-Os")
endif()
if(NOT CMAKE_CXX_FLAGS_RELEASE)
  set(CMAKE_CXX_FLAGS_RELEASE           "-O3 -DNDEBUG")
endif()
if(NOT CMAKE_CXX_FLAGS_RELWITHDEBINFO)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO    "-O2 -g -fno-omit-frame-pointer")
endif()
# set build default build type to RelWithDebInfo if not set
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# set build default build type to RelWithDebInfo if not set
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# message(STATUS "CXX_FLAGS: ${ALL_CXX_FLAGS}")
string(TOUPPER "${CMAKE_BUILD_TYPE}" UPPER_BUILD_TYPE)
set(ALL_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_BUILD_TYPE}}")

if(NOT UCR_NO_SUMMARY)
  message(STATUS
        "\n====================|  Build Summary  |===================="
        "\n"
        "\nBuild type:            ${CMAKE_BUILD_TYPE}"
        "\n"
        "\nCXX:                   ${CMAKE_CXX_COMPILER}"
        "\nCXXFLAGS:              ${ALL_CXX_FLAGS}"
        "\nLDFLAGS:               ${LDFLAGS}"
        "\n"
        "\nSource directory:      ${CMAKE_CURRENT_SOURCE_DIR}"
        "\nBuild directory:       ${CMAKE_CURRENT_BINARY_DIR}"
        "\nExecutable path:       ${EXECUTABLE_OUTPUT_PATH}"
        "\nLibrary path:          ${LIBRARY_OUTPUT_PATH}"
        "\nInstall prefix:        ${CMAKE_INSTALL_PREFIX}"
        "\nGenerator:             ${CMAKE_GENERATOR}"
        "\n"
        "\n===========================================================\n")
endif()

include_directories(. include)

file(GLOB_RECURSE HEADERS "include/*.hpp")

add_executable(ucr_ed_orig original/UCR_ED.cpp)
target_link_libraries(ucr_ed_orig ${CMAKE_DL_LIBS})

add_executable(ucr_dtw_orig original/UCR_DTW.cpp)
target_link_libraries(ucr_dtw_orig ${CMAKE_DL_LIBS})

add_executable(ucr_ed src/ucr_ed.cpp ${HEADERS})
target_link_libraries(ucr_ed ${CMAKE_DL_LIBS})

add_executable(ucr_dtw src/ucr_dtw.cpp ${HEADERS})
target_link_libraries(ucr_dtw ${CMAKE_DL_LIBS})
