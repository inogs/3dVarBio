MODULE mpi_str
  use set_knd
  use mpi

  IMPLICIT NONE

  public

  !-------------------------------------------------------!
  !     MPI variables
  !
  !     NPE : Number of Processing Elements
  !     MyId : process number  [ 0 - NPE-1 ]
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
  !     MpiWinN3n : Window for one-sided communication on grd%n3n array
  !     MpiWinN3nAd : Window for one-sided communication on grd%n3n_ad array
  !     MpiWinO2o : Window for one-sided communication on grd%o2o array
  !     MpiWinO2oAd : Window for one-sided communication on grd%o2o_ad array
  !     NextLocalRow : size of the local number of row for the process "below" MyID
  !     
  !     Var3DCommunicator : MPI Communicator (useful for the "interaction" with ogstm)
  !
  !-------------------------------------------------------!

  integer  :: NPE, MyId, MyPosI, MyPosJ
  integer  :: ProcBottom, ProcTop
  integer  :: NumProcI, NumProcJ
  integer  :: GlobalRow, GlobalCol
  integer  :: localRow, localCol
  integer  :: GlobalRowOffset, GlobalColOffset
  integer  :: MyPair
  integer  :: MpiWinChl, MpiWinChlAd
  integer  :: MpiWinN3n, MpiWinN3nAd
  integer  :: MpiWinO2o, MpiWinO2oAd
  integer  :: NextLocalRow

  integer  :: Var3DCommunicator
  integer(KIND=MPI_OFFSET_KIND) :: MyStart(3), MyCount(3)

  ! Arrays needed for alltoallv communication
  ! X dimension
  integer, allocatable, dimension(:) :: SendCountX2D, SendCountX3D, SendDisplX2D, SendDisplX3D
  integer, allocatable, dimension(:) :: RecCountX2D, RecCountX3D, RecDisplX2D, RecDisplX3D

  ! Arrays needed for the ghost cells exchange
  REAL(r8), POINTER, DIMENSION(:,:)    ::  ChlExtended
  REAL(r8), POINTER, DIMENSION(:)        ::  SendTop, RecBottom, SendBottom, RecTop

END MODULE mpi_str
