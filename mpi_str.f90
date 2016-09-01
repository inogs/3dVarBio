MODULE mpi_str
  use mpi
  
  IMPLICIT NONE
  public
  
  
  !-------------------------------------------------------!
  !     MPI vaiables
  !
  !     size : number of processes
  !     MyRank : process number  [ 0 - size-1 ]
  !     NumProcI : number of processes along i direction
  !     NumProcJ : number of processes along j direction
  !     GlobalRow : global number of rows
  !     GlobalCol : global number of columns
  !     localRow : number of row slicing in i direction
  !     localCol : number of col slicing in j direction
  !
  !-------------------------------------------------------!
  
  INTEGER  :: size, MyRank
  INTEGER  :: NumProcI, NumProcJ, NumProcIJ
  integer  :: GlobalRow, GlobalCol
  integer  :: localRow, localCol
  
  integer(KIND=MPI_OFFSET_KIND) :: MyStart(3), MyCount(3)

END MODULE mpi_str
