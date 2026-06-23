cmake_minimum_required(VERSION 3.15)
project(M3C LANGUAGES Fortran C)

# Allow selection of standard build types (Release, Debug, RelWithDebInfo, MinSizeRel)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "C compiler: ${CMAKE_C_COMPILER}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Automatically find dependencies under active Conda environments
if(DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Found CONDA_PREFIX: $ENV{CONDA_PREFIX}, adding to search paths")
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()

# Find math libraries (BLAS/LAPACK or Intel MKL)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(WITH_MKL_DEFAULT ON)
else()
  set(WITH_MKL_DEFAULT OFF)
endif()
option(WITH_MKL "Use Intel MKL if available" ${WITH_MKL_DEFAULT})

if(WITH_MKL)
  set(MKL_INTERFACE "lp64" CACHE STRING "Use 32-bit integer LP64 MKL interface" FORCE)
  find_package(MKL CONFIG QUIET)
  if(MKL_FOUND)
    message(STATUS "Found Intel MKL CONFIG package")
    set(LAPACK_LIBRARIES MKL::MKL)
  else()
    message(STATUS "Intel MKL CONFIG package not found, falling back to FindLAPACK")
    find_package(LAPACK REQUIRED)
  endif()
else()
  find_package(LAPACK REQUIRED)
endif()

# Find FFTW3
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
  if(TARGET PkgConfig::FFTW3)
    set(FFTW3_LIBRARIES PkgConfig::FFTW3)
    set(FFTW3_FOUND TRUE)
  endif()
endif()

if(NOT FFTW3_FOUND)
  find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)
  find_library(FFTW3_LIBRARY NAMES fftw3)
  if(FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
    set(FFTW3_FOUND TRUE)
    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
  endif()
endif()

if(NOT FFTW3_FOUND)
  message(FATAL_ERROR "FFTW3 library not found.")
endif()

# Centralized Fortran Module Directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Locate matching SciFT build directory
get_filename_component(BUILD_DIR_NAME "${CMAKE_BINARY_DIR}" NAME)

if(NOT DEFINED SCIFT_HOME)
  if(DEFINED ENV{SCIFT_HOME})
    set(SCIFT_HOME "$ENV{SCIFT_HOME}")
  else()
    set(SCIFT_HOME "/home/aguirre/Develop/scift")
  endif()
endif()

if(NOT DEFINED SCIFT_BUILD_DIR)
  set(SCIFT_BUILD_DIR "${SCIFT_HOME}/${BUILD_DIR_NAME}")
endif()

if(NOT EXISTS "${SCIFT_BUILD_DIR}")
  message(FATAL_ERROR "SciFT build directory not found: ${SCIFT_BUILD_DIR}. Make sure SciFT is compiled first.")
endif()

message(STATUS "Using SciFT build directory: ${SCIFT_BUILD_DIR}")

# Create imported target for SciFT
add_library(SciFT::scift STATIC IMPORTED)
set_target_properties(SciFT::scift PROPERTIES
  IMPORTED_LOCATION "${SCIFT_BUILD_DIR}/src/libscift.a"
  INTERFACE_INCLUDE_DIRECTORIES "${SCIFT_BUILD_DIR}/modules"
)

# Process subdirectories
add_subdirectory(utils/libmsym)
add_subdirectory(src)

# Copy the utils directory to the build directory at build time
add_custom_target(copy_utils ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/utils"
    "${CMAKE_CURRENT_BINARY_DIR}/utils"
  COMMENT "Copying utils directory to build tree"
)

# --- Installation and Packaging Rules ---
include(GNUInstallDirs)

# Install README.md and M3Cvars.sh to the root of the distribution
install(FILES README.md M3Cvars.sh
  DESTINATION .
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

# Install top-level files in utils/ (scripts, etc.) excluding directories
file(GLOB UTILS_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "utils/*")
foreach(FILE_PATH ${UTILS_FILES})
  if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${FILE_PATH}")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${FILE_PATH}"
      DESTINATION utils
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )
  endif()
endforeach()

# Install doc tutorial PDFs
file(GLOB TUTORIAL_FILES "doc/tutorial/tutorial-*.pdf")
if(TUTORIAL_FILES)
  install(FILES ${TUTORIAL_FILES} DESTINATION doc)
endif()

# Copy molecule.* examples from SciFT examples directory to bin/
install(DIRECTORY "${SCIFT_BUILD_DIR}/examples/"
  DESTINATION bin
  FILES_MATCHING PATTERN "molecule.*"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

# --- CPack Configuration ---
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" M3C_VERSION_RAW)
  string(STRIP "${M3C_VERSION_RAW}" M3C_VERSION)
else()
  set(M3C_VERSION "1.0.0")
endif()

set(CPACK_PACKAGE_NAME "M3C")
set(CPACK_PACKAGE_VERSION "${M3C_VERSION}")
set(CPACK_PACKAGE_VENDOR "M3C Developers")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "M3C Project")
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_FILE_NAME "M3C-v${M3C_VERSION}-${CMAKE_Fortran_COMPILER_ID}-${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY ON)

include(CPack)

