ecbuild_add_executable (TARGET cube_to_target.x SOURCES cube_to_target.F90 shr_kind_mod.F90 reconstruct.F90 remap.F90 rot.F90 shared_vars.F90 subgrid_topo_ana.F90 ridge_ana.F90)
target_link_libraries (cube_to_target.x PRIVATE MPI::MPI_Fortran NetCDF::NetCDF_Fortran)
