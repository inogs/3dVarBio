subroutine mynode()
  ! ---------------------------------------------------------------------
  ! 
  !                        routine mynode
  !                      ******************
  ! 
  !   Purpose :
  !   ---------
  !      Massively parallel processors
  !      Find processor unit
  ! 
  !    Input :
  !    -----
  !       argument                :
  ! 
  !    Modifications:
  !    --------------
  !        original  : 93-09 (M. Imbard)
  !        additions : 96-05 (j. Escobar)
  !        additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
  !                           SHMEM and MPI versions
  ! -----------------------------------------------------------------------
  
  
  use mpi_str
  
  ! -----------------------------------------------------------------------
  
  implicit none

  ! 
  !  MPI VERSION
  ! 
  !         -------------
  !         Enroll in MPI
  !         -------------
  !

  INTEGER ierr
  CALL mpi_init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, MyRank,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, size,ierr)

  !*******************************************
  !
  ! read domain decomposition files
  ! some parts of this code are copied from
  ! src/General/parini.F subroutine within ogstm package
  !
  !*******************************************
  
  call COUNTLINE ('Dom_Dec_jpi.ascii', NumProcI)
  call COUNTWORDS('Dom_Dec_jpi.ascii', NumProcJ)
  
  NumProcI = 1
  NumProcJ = Size !2
  
  MyPosI = mod(MyRank, NumProcI)
  MyPosJ = MyRank / NumProcI
  
  if(NumProcI .gt. 1 .and. NumProcJ .gt. 1) then
     ProcTop  = MyRank - 1
     if(mod(MyRank, NumProcI)-1 .lt. 0) ProcTop = MPI_PROC_NULL
     ProcBottom = MyRank + 1
     if(mod(MyRank, NumProcI)+1 .ge. NumProcI) ProcBottom = MPI_PROC_NULL
     ProcRight = MyRank + NumProcI
     if(ProcRight .ge. size) ProcRight = MPI_PROC_NULL
     ProcLeft = MyRank - NumProcI
     if(ProcLeft .lt. 0) ProcLeft = MPI_PROC_NULL
  else if(NumProcJ .gt. 1) then
     ProcLeft  = MyRank - 1
     if(ProcLeft .lt. 0) ProcLeft = MPI_PROC_NULL
     ProcRight = MyRank + 1
     if(ProcRight .ge. NumProcJ) ProcRight = MPI_PROC_NULL
     ProcBottom = MPI_PROC_NULL
     ProcTop    = MPI_PROC_NULL
  else if(NumProcI .gt. 1) then
     ProcTop  = MyRank - 1
     if(ProcTop .lt. 0) ProcTop = MPI_PROC_NULL
     ProcBottom = MyRank + 1
     if(ProcBottom .ge. NumProcI) ProcBottom = MPI_PROC_NULL
     ProcLeft  = MPI_PROC_NULL
     ProcRight = MPI_PROC_NULL
  else
     print*, ""
     print*, "You are using a single MPI Process!"
     ProcLeft   = MPI_PROC_NULL
     ProcRight  = MPI_PROC_NULL
     ProcTop    = MPI_PROC_NULL
     ProcBottom = MPI_PROC_NULL
  end if

  call MPI_Comm_split(MPI_COMM_WORLD, MyPosI, MyRank, CommSliceY, ierr)
  call MPI_Comm_split(MPI_COMM_WORLD, MyPosJ, MyRank, CommSliceX, ierr)

  ALLOCATE(SendCountX2D(NumProcI), SendCountX4D(NumProcI))
  ALLOCATE(SendDisplX2D(NumProcI), SendDisplX4D(NumProcI))
  ALLOCATE(RecCountX2D(NumProcI), RecCountX4D(NumProcI))
  ALLOCATE(RecDisplX2D(NumProcI), RecDisplX4D(NumProcI))

  ALLOCATE(SendCountY2D(NumProcJ), SendCountY4D(NumProcJ))
  ALLOCATE(SendDisplY2D(NumProcJ), SendDisplY4D(NumProcJ))
  ALLOCATE(RecCountY2D(NumProcJ), RecCountY4D(NumProcJ))
  ALLOCATE(RecDisplY2D(NumProcJ), RecDisplY4D(NumProcJ))
  
  ! write(*,*) "MyRank", MyRank, "PosI", MyPosI, "PosJ", MyPosJ, "Left", ProcLeft, "Right", ProcRight, "Top", ProcTop, "Bottom", ProcBottom

  if(NumProcI * NumProcJ .ne. size) then
     if(MyRank .eq. 0) then
        WRITE(*,*) ""
        WRITE(*,*) " Error: gridX * gridY != nproc "
        WRITE(*,*) " Exit "
        WRITE(*,*) ""
     end if
     call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  end if

  NumProcIJ = NumProcI*NumProcJ
  ! jpreci = 1
  ! jprecj = 1


  if(MyRank .eq. 0) then
     WRITE(*,*) ' '
     WRITE(*,*) 'Dom_Size'
     WRITE(*,*) ' '
     WRITE(*,*) ' number of processors following i : NumProcI   = ', NumProcI
     WRITE(*,*) ' number of processors following j : NumProcJ   = ', NumProcJ
     WRITE(*,*) ' '
     WRITE(*,*) ' local domains : < or = NumProcI x NumProcJ number of processors   = ', NumProcIJ
     ! WRITE(*,*) ' number of lines for overlap  jpreci   = ',jpreci
     ! WRITE(*,*) ' number of lines for overlap  jprecj   = ',jprecj
     WRITE(*,*) ' '
  endif

end subroutine mynode

subroutine mpi_sync

  use mpi_str
  
  implicit none

  INTEGER :: ierror

  CALL mpi_barrier(MPI_COMM_WORLD, ierror)

end subroutine mpi_sync


subroutine mpi_stop

  use mpi_str

  implicit none

  INTEGER info

  CALL mpi_finalize(info)

end subroutine mpi_stop

! **************************************************************
SUBROUTINE COUNTLINE(FILENAME,LINES)
  implicit none
  character FILENAME*(*)
  integer lines
  integer TheUnit
  
  TheUnit = 326
  
  lines=0
  OPEN(UNIT=TheUnit,file=FILENAME,status='old')
  DO WHILE (.true.)
     read(TheUnit, *, END=21)
     lines = lines+1
  ENDDO
  
21 CLOSE(TheUnit)

END SUBROUTINE COUNTLINE

! **************************************************************
SUBROUTINE COUNTWORDS(filename,n)
  IMPLICIT NONE
  CHARACTER*(*) filename
  INTEGER N
  ! local
  INTEGER I
  CHARACTER(LEN=1024) str, str_blank
  
  
  open(unit=21,file=filename, form='formatted')
  read(21,'(A)') str
  close(21)
  
  str_blank=' '//trim(str)
  N=0
  do i = 1,len(trim(str))
     if ((str_blank(i:i).eq.' ').and.(str_blank(i+1:i+1).ne.' ') )  N=N+1
  enddo
  
END SUBROUTINE COUNTWORDS
