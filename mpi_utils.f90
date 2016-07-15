SUBROUTINE mynode()
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
  
  IMPLICIT NONE

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

END SUBROUTINE mynode

SUBROUTINE mpi_sync

  USE myalloc_mpi
  
  IMPLICIT NONE

  INTEGER :: ierror

  CALL mpi_barrier(MPI_COMM_WORLD, ierror)

END SUBROUTINE mpi_sync


SUBROUTINE mpi_stop

  USE myalloc_mpi

  IMPLICIT NONE

  INTEGER info

  ! CALL mpi_sync

  CALL mpi_finalize(info)


END SUBROUTINE mpi_stop
