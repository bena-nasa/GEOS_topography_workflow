#--------------------
ecbuild_add_executable (TARGET filter_topo.x SOURCES filter_topo.F90 utils.F90)
target_link_libraries (filter_topo.x PRIVATE MPI::MPI_Fortran NetCDF::NetCDF_Fortran)
# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_link_libraries(filter_topo.x PRIVATE OpenMP::OpenMP_Fortran)
endif ()
