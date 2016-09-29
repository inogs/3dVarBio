MODULE mpi_str
  use set_knd
  use mpi
  
  IMPLICIT NONE
  public
  
  
  !-------------------------------------------------------!
  !     MPI vaiables
  !
  !     size : number of processes
  !     MyRank : process number  [ 0 - size-1 ]
  !     MyColRank : process rank in a column communicator
  !     MyRowRank : process rank in a row communicator
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
  !     GlobalRowOffset : offset needed to read grd%global_msk
  !     GlobalColOffset : offset needed to read grd%global_msk
  !
  !-------------------------------------------------------!
  
  integer  :: size, MyRank, MyPosI, MyPosJ
  integer  :: MyColRank, MyRowRank
  integer  :: ProcLeft, ProcRight, ProcBottom, ProcTop
  integer  :: NumProcI, NumProcJ, NumProcIJ
  integer  :: GlobalRow, GlobalCol
  integer  :: localRow, localCol
  integer  :: GlobalRowOffset, GlobalColOffset
  
  integer  :: CommSliceY, CommSliceX
  integer(KIND=MPI_OFFSET_KIND) :: MyStart(3), MyCount(3)

  ! Arrays needed for alltoallv communication
  ! X dimension
  integer, allocatable, dimension(:) :: SendCountX2D, SendCountX4D, SendDisplX2D, SendDisplX4D
  integer, allocatable, dimension(:) :: RecCountX2D, RecCountX4D, RecDisplX2D, RecDisplX4D
  ! Y dimension
  integer, allocatable, dimension(:) :: SendCountY2D, SendCountY4D, SendDisplY2D, SendDisplY4D
  integer, allocatable, dimension(:) :: RecCountY2D, RecCountY4D, RecDisplY2D, RecDisplY4D

  ! Arrays needed for the ghost cells exchange
  REAL(r8), POINTER, DIMENSION(:,:,:)    ::  ChlExtended
  REAL(r8), POINTER, DIMENSION(:)        ::  SendLeft, RecRight, SendRight, RecLeft
  REAL(r8), POINTER, DIMENSION(:)        ::  SendTop, RecBottom, SendBottom, RecTop

  REAL(r8), POINTER, DIMENSION(:,:,:,:)  ::  ChlExtended4D, ChlExtendedAD_4D
  REAL(r8), POINTER, DIMENSION(:,:)      ::  SendLeft2D, RecRight2D
  REAL(r8), POINTER, DIMENSION(:,:)      ::  SendRight2D, RecLeft2D
  REAL(r8), POINTER, DIMENSION(:,:)      ::  SendTop2D, RecBottom2D
  REAL(r8), POINTER, DIMENSION(:,:)      ::  SendBottom2D, RecTop2D
  

END MODULE mpi_str
