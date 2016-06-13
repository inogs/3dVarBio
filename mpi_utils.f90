SUBROUTINE mynode()
  ! INTEGER FUNCTION mynode()
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
  
  
  USE myalloc_mpi
  
  ! -----------------------------------------------------------------------
  
  IMPLICIT NONE
  
#ifdef key_mpp_mpi
  ! 
  !  MPI VERSION
  ! 
  !         -------------
  !         Enroll in MPI
  !         -------------
  !
  INTEGER ierr
  
  CALL mpi_init(ierr)
  CALL mpi_comm_rank(mpi_comm_world,rank,ierr)
  CALL mpi_comm_size(mpi_comm_world,size,ierr)
  
#  else
  INTEGER mynode
  mynode=0
  RETURN
#endif
  
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
