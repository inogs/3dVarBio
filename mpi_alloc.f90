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
  integer(KIND=MPI_OFFSET_KIND) :: MyStart(2), MyCount(2)
  
END MODULE myalloc_mpi
