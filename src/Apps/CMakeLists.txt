#--------------------
# Copy include files that are used by other libraries.
# We could leave these in the source directory, and just broaden the search path
# in the other libaries, but this make it explicit which aspects are externally
# used.

ecbuild_add_executable (TARGET generate_scrip_cube.x SOURCES generate_scrip_cube.F90 geompack.F90)
target_link_libraries (generate_scrip_cube.x PRIVATE MAPL MPI::MPI_Fortran esmf)
target_include_directories (generate_scrip_cube.x PRIVATE $<BUILD_INTERFACE:${MAPL_SOURCE_DIR}/include>)
# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_link_libraries(generate_scrip_cube.x PRIVATE OpenMP::OpenMP_Fortran)
endif ()

