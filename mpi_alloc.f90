MODULE myalloc_mpi
  use mpi
  
  IMPLICIT NONE
  public
  
  
!#include <mpif.h>
  
  !-----------------------------------------------!
  !     MPI vaiables
  !
  !     size : number of process
  !     rank : process number  [ 0 - size-1 ]
  !-----------------------------------------------!
  
  INTEGER  :: size,rank
  ! INTEGER  :: comm ! MPI_COMM
  
END MODULE myalloc_mpi
