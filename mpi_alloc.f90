MODULE myalloc_mpi
  use mpi
  
  IMPLICIT NONE
  public
  
  
  !-----------------------------------------------!
  !     MPI vaiables
  !
  !     size : number of process
  !     rank : process number  [ 0 - size-1 ]
  !-----------------------------------------------!
  
  INTEGER  :: size,rank
  
END MODULE myalloc_mpi
