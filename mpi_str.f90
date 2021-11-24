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
  integer, allocatable, dimension(:) :: SendCountX3D_chl, SendDisplX3D_chl
  integer, allocatable, dimension(:) :: RecCountX2D, RecCountX3D, RecDisplX2D, RecDisplX3D
  integer, allocatable, dimension(:) :: RecCountX3D_chl, RecDisplX3D_chl

  ! Arrays needed for the ghost cells exchange
  REAL(r8), POINTER, DIMENSION(:,:)    ::  ChlExtended
  REAL(r8), POINTER, DIMENSION(:)        ::  SendTop, RecBottom, SendBottom, RecTop
  REAL(r8), ALLOCATABLE, DIMENSION(:,:)      ::  SendTop_2d, RecBottom_2d
  REAL(r8), ALLOCATABLE, DIMENSION(:,:)      ::  SendBottom_2d, RecTop_2d
  REAL(r8), ALLOCATABLE,  DIMENSION(:,:,:)    ::  ChlExtended_3d, N3nExtended_3d, O2oExtended_3d


CONTAINS

  SUBROUTINE EXTEND_2D(INPUT, my_km, OUTPUT_Extended)
  use set_knd
  use grd_str
  use obs_str
  use drv_str
  use bio_str
  IMPLICIT NONE
  INTEGER(i4)   ::  my_km
  REAL(r8), DIMENSION(grd%im  , grd%jm  , my_km), INTENT(IN)  ::  INPUT
  REAL(r8), DIMENSION(grd%im+1, grd%jm  , my_km), INTENT(OUT) ::  OUTPUT_Extended


  INTEGER(i4)   ::  i, j, k, kk
  INTEGER   :: MyTag
  INTEGER   :: ReqTop, ReqBottom, ierr
  INTEGER   :: StatBottom(MPI_STATUS_SIZE)


      ! Filling array to send
      do k=1,my_km
       do j=1,grd%jm
         SendTop_2d(j,k)  = INPUT(1,j,k)
       end do
      end do

      do k=my_km+1,grd%km
       SendTop_2d(:,k)  = 0
      end do

      MyTag = 42
      RecBottom_2d(:,:) = 0

      call MPI_Isend(SendTop_2d, grd%jm*grd%km, MPI_REAL8, ProcTop, MyTag, &
           Var3DCommunicator, ReqTop, ierr)
      call MPI_Irecv(RecBottom_2d, grd%jm*grd%km, MPI_REAL8, ProcBottom, MyTag, &
           Var3DCommunicator, ReqBottom, ierr)

      do k=1,my_km
          do j=1,grd%jm
             do i=1,grd%im
                OUTPUT_Extended(i,j,k) = INPUT(i,j,k)
             end do
          end do
      end do


      call MPI_Wait(ReqBottom, StatBottom, ierr)
      do k=1,my_km
        do j=1,grd%jm
          OUTPUT_Extended(grd%im+1,j,k) = RecBottom_2d(j,k)
        end do
      end do



  END SUBROUTINE EXTEND_2D


  SUBROUTINE ADD_PREVCORE_CONTRIB(INPUT_Extended,my_km,OUTPUT,INIT_2d)
  ! ADDS PREVIOUS CORE CONTRIBUTION to the sum of an _ad variable
  !  used in obs_arg_ad
  !
  !  THEORY

  !  OUTPUT = INIT + contr(i) + contr(i+1)
  !  but without ghost cell we have
  ! OUTPUT(1)    = INIT + contr(1)
  ! OUTPUT(im+1) = INIT + contr(im+1)
  ! i.e one single contribution

  ! im+1 is position 1 of following core,
  ! then we can receive contr(im+1) in RecTop_2d

  ! In order to have  OUTPUT = INIT + contr(i) + contr(i+1) we'll do
  ! OUTPUT(1) = INIT + contrib(1) + INIT + contr(i+1) - INIT
  !             OUTPUT(1)         + Rec_top           - INIT

  use set_knd
  use grd_str
  use obs_str
  use drv_str
  use bio_str
  IMPLICIT NONE
  INTEGER(i4)   ::  my_km
  REAL(r8), DIMENSION(grd%im+1, grd%jm  , my_km), INTENT(IN)   ::  INPUT_Extended
  REAL(r8), DIMENSION(grd%im  , grd%jm  , my_km), INTENT(OUT)  ::  OUTPUT
  REAL(r8), DIMENSION(          grd%jm  , my_km), INTENT(IN)   ::  INIT_2d

  INTEGER   :: ReqBottom, ReqTop, ierr
  INTEGER(i4)   ::  i, j, k, kk
  INTEGER   :: MyTag
  INTEGER   :: StatTop(MPI_STATUS_SIZE)

  do k=1,my_km
    do j=1,grd%jm
      SendBottom_2d(j,k) = INPUT_Extended(grd%im+1,j,k)
    end do
  end do


  do k=my_km+1,grd%km
    SendBottom_2d(:,k) = 0
  end do

  MyTag = 42
  RecTop_2d(:,:) = 0

  call MPI_Isend(SendBottom_2d, grd%jm*grd%km, MPI_REAL8, ProcBottom, MyTag, &
       Var3DCommunicator, ReqBottom, ierr)
  call MPI_Irecv(RecTop_2d, grd%jm*grd%km, MPI_REAL8, ProcTop, MyTag, &
       Var3DCommunicator, ReqTop, ierr)

   OUTPUT = INPUT_Extended(1:grd%im,:,:)

  call MPI_Wait(ReqTop, StatTop, ierr)

  do k=1,my_km
    do j=1,grd%jm
      OUTPUT(1,j,k) = OUTPUT(1,j,k) + RecTop_2d(j,k) - INIT_2d(j,k)
    end do
  end do

  END SUBROUTINE ADD_PREVCORE_CONTRIB

END MODULE mpi_str
