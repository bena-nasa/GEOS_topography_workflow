ecbuild_add_executable (TARGET bin_to_cube.x SOURCES bin_to_cube.F90 shr_kind_mod.F90)
target_link_libraries (bin_to_cube.x PRIVATE MPI::MPI_Fortran NetCDF::NetCDF_Fortran)
install(PROGRAMS  landm_coslat.nc DESTINATION bin)
