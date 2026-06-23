cmake_minimum_required(VERSION 3.5)

# Set C99 standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Compile flags and include path
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${CMAKE_CURRENT_SOURCE_DIR} -static -static-libgcc")

add_library(msym STATIC
  basis_function.c
  character_table.c
  context.c
  debug.c
  elements.c
  equivalence_set.c
  geometry.c
  linalg.c
  msym.c
  msym_error.c
  permutation.c
  point_group.c
  rsh.c
  subspace.c
  symmetrize.c
  symmetry.c
  symop.c
)

add_executable(molecule.symmetrize molecule.symmetrize.c)
target_link_libraries(molecule.symmetrize msym m)

set_target_properties(molecule.symmetrize PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/utils"
)

install(TARGETS molecule.symmetrize DESTINATION utils)