program ocean_var
use filenames
implicit none

call SETFILENAMES
call my_mpi_init
call oceanvar

end program ocean_var
