cmake_minimum_required(VERSION 2.8)
project(msym)
#add_definitions(-std=c99 -I. -static -static-libstdc++ -static-libgcc)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -I. -static -static-libgcc -static-libstdc++")

add_library(msym STATIC
basis_function.c   character_table.h  debug.c     elements.h         geometry.h  msym.c        msym.h         point_group.c  rsh.h       symmetrize.c  symmetry.h
basis_function.h   context.c          debug.h     equivalence_set.c  linalg.c    msym_error.c  permutation.c  point_group.h  subspace.c  symmetrize.h  symop.c
character_table.c  context.h          elements.c  equivalence_set.h  geometry.c  linalg.h    msym_error.h  permutation.h  rsh.c          subspace.h  symmetry.c    symop.h
)

add_executable(molecule.symmetrize molecule.symmetrize.c)
target_link_libraries(molecule.symmetrize msym)
target_link_libraries(molecule.symmetrize m)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../)
install(TARGETS molecule.symmetrize DESTINATION bin)