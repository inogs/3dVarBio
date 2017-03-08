program ocean_var
use filenames
use mpi_str
implicit none

call SETFILENAMES

! init of the MPI environment
! in case of standalone usage
call my_mpi_init

call oceanvar

! finalizing the MPI environment
call MPI_Comm_free(Var3DCommunicator, ierr)
call mpi_stop

end program ocean_var
