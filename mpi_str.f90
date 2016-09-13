MODULE mpi_str
  use mpi
  
  IMPLICIT NONE
  public
  
  
  !-------------------------------------------------------!
  !     MPI vaiables
  !
  !     size : number of processes
  !     MyRank : process number  [ 0 - size-1 ]
  !     MyPosI : rank on communicator along row direction
  !     MyPosJ : rank on communicator along column direction
  !     ProcLeft : rank of process to my left side
  !     ProcRight: rank of process to my right side
  !     ProcBottom: rank of process under me
  !     ProcTop: rank of process on top of me
  !     NumProcI : number of processes along i direction
  !     NumProcJ : number of processes along j direction
  !     GlobalRow : global number of rows
  !     GlobalCol : global number of columns
  !     localRow : number of row slicing in i direction
  !     localCol : number of col slicing in j direction
  !
  !-------------------------------------------------------!
  
  integer  :: size, MyRank, MyPosI, MyPosJ
  integer  :: ProcLeft, ProcRight, ProcBottom, ProcTop
  integer  :: NumProcI, NumProcJ, NumProcIJ
  integer  :: GlobalRow, GlobalCol
  integer  :: localRow, localCol
  
  integer(KIND=MPI_OFFSET_KIND) :: MyStart(3), MyCount(3)

END MODULE mpi_str
