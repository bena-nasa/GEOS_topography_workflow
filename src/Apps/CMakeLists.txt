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

ecbuild_add_executable (TARGET convert_bin_to_netcdf.x SOURCES convert_bin_to_netcdf.F90)
target_link_libraries (convert_bin_to_netcdf.x PRIVATE MAPL MPI::MPI_Fortran)
#target_include_directories (convert_bin_to_netcdf.x PRIVATE $<BUILD_INTERFACE:${MAPL_SOURCE_DIR}/include>)
# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_link_libraries(convert_bin_to_netcdf.x PRIVATE OpenMP::OpenMP_Fortran)
endif ()

ecbuild_add_executable (TARGET convert_to_gmao_output.x SOURCES convert_to_gmao_output.F90)
target_link_libraries (convert_to_gmao_output.x PRIVATE MAPL MPI::MPI_Fortran)
#target_include_directories (convert_to_gmao_output.x PRIVATE $<BUILD_INTERFACE:${MAPL_SOURCE_DIR}/include>)
# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_link_libraries(convert_to_gmao_output.x PRIVATE OpenMP::OpenMP_Fortran)
endif ()

file(COPY  mosaic_grid.py DESTINATION bin)
file(COPY  mosaic_topo.py DESTINATION bin)
file(COPY  combineMosaic.py DESTINATION bin)

esma_add_subdirectories(
   bin_to_cube
   cube_to_target
   )
