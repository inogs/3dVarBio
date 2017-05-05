program ocean_var
use filenames
implicit none

integer :: ierr

call SETFILENAMES

! init of the MPI environment
! in case of standalone usage
call var3d_mpi_init

call oceanvar

! finalizing the MPI environment
call mpi_stop

end program ocean_var
