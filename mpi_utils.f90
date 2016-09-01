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
  NumProcJ = 2

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
