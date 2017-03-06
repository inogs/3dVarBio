subroutine my_mpi_init()

  use mpi_str
  implicit none

  integer :: ierr, zero

  CALL mpi_init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, MyRank,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, size,ierr)

  zero = 0
  call MPI_Comm_split(MPI_COMM_WORLD, zero, MyRank, MyCommWorld, ierr)

end subroutine my_mpi_init

subroutine my_3dvar_node()
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

  NumProcI = Size
  NumProcJ = 1

  MyPosI = mod(MyRank, NumProcI)
  MyPosJ = MyRank / NumProcI

  if(NumProcI .gt. 1) then
     ProcTop  = MyRank - 1
     if(ProcTop .lt. 0) ProcTop = MPI_PROC_NULL
     ProcBottom = MyRank + 1
     if(ProcBottom .ge. NumProcI) ProcBottom = MPI_PROC_NULL
  else
     print*, ""
     print*, "You are using a single MPI Process!"
     ProcTop    = MPI_PROC_NULL
     ProcBottom = MPI_PROC_NULL
  end if

  call MPI_TYPE_CONTIGUOUS(2, MPI_REAL8, MyPair, ierr)
  call MPI_TYPE_COMMIT(MyPair, ierr)

  ALLOCATE(SendCountX2D(NumProcI), SendCountX4D(NumProcI))
  ALLOCATE(SendDisplX2D(NumProcI), SendDisplX4D(NumProcI))
  ALLOCATE(RecCountX2D(NumProcI), RecCountX4D(NumProcI))
  ALLOCATE(RecDisplX2D(NumProcI), RecDisplX4D(NumProcI))

  ! print for debug 
  ! write(*,*) "MyRank", MyRank, "PosI", MyPosI, "PosJ", MyPosJ, "Left", ProcLeft, "Right", ProcRight, "Top", ProcTop, "Bottom", ProcBottom

  if(NumProcI * NumProcJ .ne. size) then
     if(MyRank .eq. 0) then
        WRITE(*,*) ""
        WRITE(*,*) " Error: gridX * gridY != nproc "
        WRITE(*,*) " Exit "
        WRITE(*,*) ""
     end if
     call MPI_Abort(MyCommWorld, -1, ierr)
  end if

  if(MyRank .eq. 0) then
     WRITE(*,*) ' '
     WRITE(*,*) 'Dom_Size'
     WRITE(*,*) ' '
     WRITE(*,*) ' number of processors following i : NumProcI   = ', NumProcI
     WRITE(*,*) ' number of processors following j : NumProcJ   = ', NumProcJ
     WRITE(*,*) ' '
  endif

end subroutine my_3dvar_node

subroutine mpi_sync

  use mpi_str

  implicit none

  INTEGER :: ierror

  CALL mpi_barrier(MyCommWorld, ierror)

end subroutine mpi_sync


subroutine mpi_stop

  use mpi_str

  implicit none

  INTEGER info

  CALL mpi_finalize(info)

end subroutine mpi_stop
