program ocean_var
use filenames
implicit none

call SETFILENAMES

! init of the MPI environment
! in case of standalone usage
call my_mpi_init

call oceanvar

! finalizing the MPI environment
call mpi_stop

end program ocean_var
