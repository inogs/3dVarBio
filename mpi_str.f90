MODULE mpi_str
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
  integer  :: jpiglo, jpjglo
  integer(KIND=MPI_OFFSET_KIND) :: MyStart(3), MyCount(3)
  
END MODULE mpi_str
