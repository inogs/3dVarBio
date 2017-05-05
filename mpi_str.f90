MODULE mpi_str
  use set_knd
  use mpi

  IMPLICIT NONE
  !include "mpif.h"
  public

  !-------------------------------------------------------!
  !     MPI vaiables
  !
  !     NPE : Number of Processing Elements
  !     MyId : process number  [ 0 - NPE-1 ]
  !     MyColRank : process rank in a column communicator
  !     MyRowRank : process rank in a row communicator
  !     MyPosI : rank on communicator along row direction
  !     MyPosJ : rank on communicator along column direction
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
  !     MpiWinChl : Window for one-sided communication on grd%chl array
  !     MpiWinChlAd : Window for one-sided communication on grd%chl_ad array
  !     NextLocalRow : NPE of the local number of row for the process "below" me
  !
  !-------------------------------------------------------!

  integer  :: NPE, MyId, MyPosI, MyPosJ
  integer  :: MyColRank, MyRowRank
  integer  :: ProcBottom, ProcTop
  integer  :: NumProcI, NumProcJ
  integer  :: GlobalRow, GlobalCol
  integer  :: localRow, localCol
  integer  :: GlobalRowOffset, GlobalColOffset
  integer  :: MyPair
  integer  :: MpiWinChl, MpiWinChlAd
  integer  :: NextLocalRow

  integer  :: Var3DCommunicator
  integer(KIND=MPI_OFFSET_KIND) :: MyStart(3), MyCount(3)

  ! Arrays needed for alltoallv communication
  ! X dimension
  integer, allocatable, dimension(:) :: SendCountX2D, SendCountX4D, SendDisplX2D, SendDisplX4D
  integer, allocatable, dimension(:) :: RecCountX2D, RecCountX4D, RecDisplX2D, RecDisplX4D

  ! Arrays needed for the ghost cells exchange
  REAL(r8), POINTER, DIMENSION(:,:,:)    ::  ChlExtended
  REAL(r8), POINTER, DIMENSION(:)        ::  SendTop, RecBottom, SendBottom, RecTop

END MODULE mpi_str
