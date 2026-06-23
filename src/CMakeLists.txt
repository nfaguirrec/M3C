set(M3C_lib_SOURCES
  AzizSlamanPotential.f90
  Fragment.f90
  FragmentsDB.f90
  FragmentsListBase.f90
  FragmentsList.f90
  GOptionsM3C.f90
  MarkovChain.f90
  ModelPotential.f90
  NNLS.f90
  Reactor.f90
)

# Build the main M3C library (static)
add_library(M3C_lib STATIC ${M3C_lib_SOURCES})
set_target_properties(M3C_lib PROPERTIES 
  OUTPUT_NAME M3C
  POSITION_INDEPENDENT_CODE ON
)

# Compiler compiler options
target_compile_options(M3C_lib PUBLIC
  $<$<Fortran_COMPILER_ID:GNU>:-w -cpp -fdec-structure -ffree-line-length-none>
  $<$<Fortran_COMPILER_ID:Intel>:-w -fpp>
  $<$<Fortran_COMPILER_ID:IntelLLVM>:-w -fpp>
)

# Include directories
target_include_directories(M3C_lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link libraries
target_link_libraries(M3C_lib PUBLIC
  SciFT::scift
  ${LAPACK_LIBRARIES}
  ${FFTW3_LIBRARIES}
)

# Compile executables
add_executable(M3C M3C.f90)
target_link_libraries(M3C PRIVATE M3C_lib)

add_executable(M3CBR M3CBR.f90)
target_link_libraries(M3CBR PRIVATE M3C_lib)

add_executable(M3CfitBR M3CfitBR.f90)
target_link_libraries(M3CfitBR PRIVATE M3C_lib)

add_executable(molecule.ZPE molecule.ZPE.f90)
target_link_libraries(molecule.ZPE PRIVATE M3C_lib)

# Install main executables
install(TARGETS M3C M3CBR M3CfitBR molecule.ZPE
  RUNTIME DESTINATION bin
)

