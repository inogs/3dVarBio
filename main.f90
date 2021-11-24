program ocean_var
use filenames
use da_params
implicit none

integer :: ierr

call SETFILENAMES
call SET_DA_PARAMS

! init of the MPI environment
! in case of standalone usage
call var3d_mpi_init

call oceanvar

! finalizing the MPI environment
call clean_da_params
call mpi_stop

end program ocean_var
