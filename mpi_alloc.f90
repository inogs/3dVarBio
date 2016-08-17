MODULE myalloc_mpi
  use mpi
  
  IMPLICIT NONE
  public
  
  
  !-----------------------------------------------!
  !     MPI vaiables
  !
  !     size : number of process
  !     MyRank : process number  [ 0 - size-1 ]
  !-----------------------------------------------!
  
  INTEGER  :: size, MyRank
  INTEGER  :: jpni, jpnj, jpnij
  
END MODULE myalloc_mpi
