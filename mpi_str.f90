MODULE mpi_str
  use set_knd
  use mpi

  IMPLICIT NONE
  public


  !-----------------------------------------------------------------!
  !     MPI variables
  !
  !     Size : number of processes
  !     MyRank : process number  [ 0 - size-1 ]
  !     LevSize : size of block levels to send to slaves
  !     NLevels : number of blocks (Total number of levels / LevSize)
  !     LevRest : relative rest (Total number of levels % Size)
  !
  !-----------------------------------------------------------------!

  integer  :: Size, MyRank
  integer  :: LevSize, NLevels, LevRest


END MODULE mpi_str
