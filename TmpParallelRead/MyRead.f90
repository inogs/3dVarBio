program MyPnetCDF
  
  use mpi
  use pnetcdf

  implicit none
  
  integer :: ierr, MyID, size, cmode, ncid
  integer :: jpni, jpnj, jpnij, jpreci, jprecj, jpkb
  integer :: jpiglo, jpjglo, jpk, DimId, VarId
  real(4), allocatable :: values(:,:)

  integer, allocatable :: ilcit(:,:), ilcjt(:,:)
  integer :: ji, jj, jpi, jpj, nn !, jpij, jpim1, jpjm1, jpkm1, jpkbm1
  integer :: MyRest, RealOffset, dimd, xtype, ndims
  integer, allocatable :: dimids(:)
  
  integer(KIND=MPI_OFFSET_KIND) MyOffset
  integer(KIND=MPI_OFFSET_KIND) MyStart(2), MyCount(2), MyTest(2)
  character(LEN=256) :: filename, MyVar

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyID, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

  ! init check
  write(*,*) "Hello World by process ", myid

  ! read domain decomposition files
  call COUNTLINE ('Dom_Dec_jpi.ascii', jpni)
  call COUNTWORDS('Dom_Dec_jpi.ascii', jpnj)
  
  jpnij = jpni*jpnj
  jpreci = 1
  jprecj = 1
  


  if(MyID .eq. 0) then
     WRITE(*,*) ' '
     WRITE(*,*) 'Dom_Size'
     WRITE(*,*) ' '
     WRITE(*,*) ' number of processors following i : jpni   = ', jpni
     WRITE(*,*) ' number of processors following j : jpnj   = ', jpnj
     WRITE(*,*) ' '
     WRITE(*,*) ' local domains : < or = jpni x jpnj number of processors   = ', jpnij
     WRITE(*,*) ' number of lines for overlap  jpreci   = ',jpreci
     WRITE(*,*) ' number of lines for overlap  jprecj   = ',jprecj
     WRITE(*,*) ' '
  endif
  
  filename = "meshmask_872.nc"
  cmode = NF90_NOWRITE
  ierr = nf90mpi_open(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, ncid)

  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_open', ierr)

  call MyGetDimension(ncid, 'x', MyOffset)
  jpiglo = MyOffset

  call MyGetDimension(ncid, 'y', MyOffset)
  jpjglo = MyOffset

  call MyGetDimension(ncid, 'z', MyOffset)
  jpk = MyOffset
  jpkb = jpk

  if(MyID .eq. 0) then
     WRITE(*,*) 'Dimension_Med_Grid'
     WRITE(*,*) ' '
     WRITE(*,*) ' jpiglo  : first  dimension of global domain --> i ',jpiglo
     WRITE(*,*) ' jpjglo  : second dimension of global domain --> j ',jpjglo
     WRITE(*,*) ' jpk     : number of levels           > or = jpk   ',jpk
     ! WRITE(*,*) ' jpkb    : first vertical layers where biology is active > or = jpkb   ',jpkb
     WRITE(*,*) ' WorkLoad: jpiglo / size                           ',jpiglo/ size
     WRITE(*,*) ' Rest    : mod(jpiglo, size)                       ',mod(jpiglo, size)
     WRITE(*,*) ' '
  endif

  !*******************************************

  allocate(ilcit(jpni, jpnj)) ; ilcit = huge(ilcit(1,1))
  allocate(ilcjt(jpni, jpnj)) ; ilcjt = huge(ilcjt(1,1))
  
  open(3333,file='Dom_Dec_jpi.ascii', form='formatted')
  open(3334,file='Dom_Dec_jpj.ascii', form='formatted')
  
  read(3333,*) ((ilcit(ji,jj), jj=1,jpnj),ji=1,jpni)
  read(3334,*) ((ilcjt(ji,jj), jj=1,jpnj),ji=1,jpni)
  
  close(3333)
  close(3334)

  do nn =1, jpni*jpnj
     if(MyId+1 .EQ. nn) then
        ji = 1 + mod(nn -1, jpni)
        jj = 1 + (nn -1)/jpni
        jpi =  ilcit(ji,jj) 
        jpj =  ilcjt(ji,jj)
     endif
  enddo
  
  ! "global" value not needed at this moment
  ! jpim1=jpi-1
  ! jpjm1=jpj-1
  ! jpkm1=jpk-1
  ! jpij=jpi*jpj
  ! jpkbm1=jpkb-1
  
  MyRest = mod(jpiglo, size)
  MyCount(2) = jpiglo / size
  RealOffset = 0
  if (MyId .lt. MyRest) then
     MyCount(2) = MyCount(2) + 1
  else
     RealOffset = MyRest
  end if

  MyStart(2) = MyCount(2) * MyID + RealOffset + 1
  MyCount(2) = MyCount(2) - 1

  MyStart(1) = 1
  MyCount(1) = jpjglo -1

  if(MyID .eq. 0) then
     write(*,*) "MyID = ", MyId, " MyStart = ", MyStart, " MyCount = ", &
          MyCount, " Sum = ", MyCount + MyStart
  else
     write(*,*) "MyID = ", MyId, " MyStart = ", MyStart, " MyCount = ", &
          MyCount, " Sum = ", MyCount +MyStart
  end if
  
  allocate(values(jpjglo, MyCount(2)+1))
  ! allocate(values(jpiglo, jpiglo))

  ! write(*,*) "check dimensions: MyID = ", MyID, " shape(values) = ", shape(values)

  ierr = nf90mpi_inq_varid(ncid, "nav_lon", VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)

  ! MyVar = "nav_lon"
  ! allocate(dimids(2))
  ! ierr = nf90mpi_inquire_variable(ncid, VarId, MyVar, xtype, ndims, dimids)
  ! if(MyID .eq.0) then
  !    write(*,*) "MyVar ", trim(MyVar), " xtype = ", xtype, " ndims = ", ndims, " dimids = ", dimids
  ! end if

  ! ierr = nfmpi_get_vara_double_all(ncid, VarId, MyStart, MyCount, values)
  ierr = nfmpi_get_vara_real_all(ncid, VarId, MyStart, MyCount, values)
  if (ierr .ne. NF90_NOERR ) call new_handle_err('nf90mpi_get_vara_real', ierr, MyId)

  ! write(*,*) "MyID = ", MyID, " shape(values) = ", shape(values)
  ! write(*,*) "MyID = ", MyID, " values = ", values

  !*******************************************

  ierr = nf90mpi_close(ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_close', ierr)

  call MPI_Finalize(ierr)

end program MyPnetCDF

subroutine MyGetDimension(ncid, name, n)
  use pnetcdf
  use mpi
  implicit none

  character name*(*)
  integer :: ncid, ierr
  integer(KIND=MPI_OFFSET_KIND) :: n
  ! integer(8) :: n
  integer dimid

  ierr = nf90mpi_inq_dimid(ncid, name, DimId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_dimid', ierr)
  ierr = nfmpi_inq_dimlen(ncid, DimId, n)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_inq_dimlen', ierr)

end subroutine MyGetDimension


subroutine handle_err(err_msg, errcode)

  use mpi
  use pnetcdf

  implicit none

  character*(*), intent(in) :: err_msg
  integer,       intent(in) :: errcode

  !local variables
  integer err

  write(*,*) 'Error: ', trim(err_msg), ' ', nf90mpi_strerror(errcode)
  call MPI_Abort(MPI_COMM_WORLD, -1, err)
  return
end subroutine handle_err

subroutine new_handle_err(err_msg, errcode, me)

  use mpi
  use pnetcdf

  implicit none

  character*(*), intent(in) :: err_msg
  integer,       intent(in) :: errcode, me

  !local variables
  integer err

  write(*,*) 'MyID = ', me, ' Error: ', trim(err_msg), ' ', nf90mpi_strerror(errcode)
  call MPI_Abort(MPI_COMM_WORLD, -1, err)
  return
end subroutine new_handle_err


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
