subroutine mynode()
  ! ---------------------------------------------------------------------
  ! 
  !                        routine mynode
  !                      ******************
  ! 
  !   Purpose :
  !   ---------
  !      Massively parallel processors
  !      Find processor unit
  ! 
  !    Input :
  !    -----
  !       argument                :
  ! 
  !    Modifications:
  !    --------------
  !        original  : 93-09 (M. Imbard)
  !        additions : 96-05 (j. Escobar)
  !        additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
  !                           SHMEM and MPI versions
  ! -----------------------------------------------------------------------
  
  
  use myalloc_mpi
  
  ! -----------------------------------------------------------------------
  
  implicit none

  ! 
  !  MPI VERSION
  ! 
  !         -------------
  !         Enroll in MPI
  !         -------------
  !

  INTEGER ierr
  CALL mpi_init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, rank,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, size,ierr)

end subroutine mynode

subroutine mpi_sync

  use myalloc_mpi
  
  implicit none

  INTEGER :: ierror

  CALL mpi_barrier(MPI_COMM_WORLD, ierror)

end subroutine mpi_sync


subroutine mpi_stop

  use myalloc_mpi

  implicit none

  INTEGER info

  ! CALL mpi_sync

  CALL mpi_finalize(info)


end subroutine mpi_stop
