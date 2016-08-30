subroutine parallel_rdgrd  
  
  use set_knd
  use drv_str
  use grd_str
  use filenames

  use mpi
  use mpi_str
  use pnetcdf
  
  implicit none

  integer(i8) :: ierr, ncid
  integer(i8) :: jpreci, jprecj
  integer(i8) :: jpiglo, jpjglo
  integer(i8) :: TmpInt, VarId
  real(r4), ALLOCATABLE          :: x3(:,:,:), x2(:,:), x1(:)

  integer, allocatable :: ilcit(:,:), ilcjt(:,:)
  integer(i8) :: ji, jj, jpi, jpj, nn
  integer(i8) :: MyRestRow, MyRestCol, OffsetRow, OffsetCol
  
  integer(KIND=MPI_OFFSET_KIND) MyOffset
  
  !
  ! open grid1.nc in read-only mode
  ierr = nf90mpi_open(MPI_COMM_WORLD, GRID_FILE, NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_open', ierr)

  !
  ! get grid dimensions
  !
  call MyGetDimension(ncid, 'im', MyOffset)
  grd%im = MyOffset
  jpiglo = MyOffset

  call MyGetDimension(ncid, 'jm', MyOffset)
  grd%jm = MyOffset
  jpjglo = MyOffset

  call MyGetDimension(ncid, 'km', MyOffset)
  grd%km = MyOffset

  if(MyRank .eq. 0) then
     write(drv%dia,*)'Grid dimensions are: ',grd%im,grd%jm,grd%km

     WRITE(*,*) 'Dimension_Med_Grid'
     WRITE(*,*) ' '
     WRITE(*,*) ' jpiglo  : first  dimension of global domain --> i ',jpiglo
     WRITE(*,*) ' jpjglo  : second dimension of global domain --> j ',jpjglo
     WRITE(*,*) ' '
  endif

  allocate(ilcit(jpni, jpnj)) ; ilcit = huge(ilcit(1,1))
  allocate(ilcjt(jpni, jpnj)) ; ilcjt = huge(ilcjt(1,1))
  
  open(3333,file='Dom_Dec_jpi.ascii', form='formatted')
  open(3334,file='Dom_Dec_jpj.ascii', form='formatted')
  
  read(3333,*) ((ilcit(ji,jj), jj=1,jpnj),ji=1,jpni)
  read(3334,*) ((ilcjt(ji,jj), jj=1,jpnj),ji=1,jpni)
  
  close(3333)
  close(3334)

  do nn =1, jpni*jpnj
     if(MyRank+1 .EQ. nn) then
        ji = 1 + mod(nn -1, jpni)
        jj = 1 + (nn -1)/jpni
        jpi =  ilcit(ji,jj) 
        jpj =  ilcjt(ji,jj)
     endif
  enddo  
  
  !*******************************************
  !
  ! PDICERBO version of the domain decomposition:
  ! the domain is divided among the processes into slices
  ! of size (jpiglo / jpni, jpjglo / jpnj).
  ! Clearly, the division is done tacking into account 
  ! rests. The only condition we need is that jpni*jpnj = NPROC
  !
  ! WARNING!!! netcdf stores data in ROW MAJOR order
  ! while here we are reading in column major order.
  ! We have to take into account this simply swapping ("ideally")
  ! the entries of MyStart and MyCount
  !
  !*******************************************
  
  MyRestCol = mod(jpiglo, jpni)
  MyRestRow = mod(jpjglo, jpnj)
  
  ! computing rests for X direction
  MyCount(1) = jpiglo / jpni
  MyCount(2) = jpjglo / jpnj
  OffsetCol = 0
  if (mod(MyRank, jpni) .lt. MyRestCol) then
     MyCount(1) = MyCount(1) + 1
     OffsetCol = mod(MyRank, jpni)
  else
     OffsetCol = MyRestCol
  end if
  
  ! computing rests for Y direction
  OffsetRow = 0
  TmpInt = MyRank / jpni
  if (TmpInt .lt. MyRestRow) then
     MyCount(2) = MyCount(2) + 1
  else
     OffsetRow = MyRestRow
  end if
  
  TmpInt = jpiglo / jpni
  MyStart(1) = TmpInt * mod(MyRank, jpni) + OffsetCol + 1
  MyCount(1) = MyCount(1)
  
  TmpInt = MyRank / jpni
  MyStart(2) = mod(MyCount(2) * TmpInt + OffsetRow, jpjglo) + 1
  MyCount(2) = MyCount(2)

  ! taking all values along k direction
  MyStart(3) = 1
  MyCount(3) = grd%km
  
  if(MyRank .eq. 0) then
     write(*,*) "MyRank = ", MyRank, " MyStart = ", MyStart, " MyCount = ", &
          MyCount, " Sum = ", MyCount + MyStart
  else
     write(*,*) "MyRank = ", MyRank, " MyStart = ", MyStart, " MyCount = ", &
          MyCount, " Sum = ", MyCount +MyStart
  end if

  grd%im = MyCount(1)
  grd%jm = MyCount(2)
  
  ! *****************************************************************************************
  ! *****************************************************************************************
  ! (almost) copy-paste from rdgrds.f90
  ! Allocate grid arrays

  ALLOCATE ( grd%reg(grd%im,grd%jm))        ; grd%reg = huge(grd%reg(1,1))
  ALLOCATE ( grd%msk(grd%im,grd%jm,grd%km)) ; grd%msk = huge(grd%msk(1,1,1))
  ALLOCATE ( grd%dep(grd%km))        ; grd%dep = huge(grd%dep(1))
  ALLOCATE ( grd%dx(grd%im,grd%jm))  ; grd%dx  = huge(grd%dx(1,1))
  ALLOCATE ( grd%dy(grd%im,grd%jm))  ; grd%dy  = huge(grd%dy(1,1))
  
  ALLOCATE ( grd%alx(grd%im,grd%jm) )         ; grd%alx  = huge(grd%alx(1,1))
  ALLOCATE ( grd%aly(grd%im,grd%jm) )         ; grd%aly  = huge(grd%aly(1,1))
  ALLOCATE ( grd%btx(grd%im,grd%jm) )         ; grd%btx  = huge(grd%btx(1,1))
  ALLOCATE ( grd%bty(grd%im,grd%jm) )         ; grd%bty  = huge(grd%bty(1,1))
  ALLOCATE ( grd%scx(grd%im,grd%jm) )         ; grd%scx  = huge(grd%scx(1,1))
  ALLOCATE ( grd%scy(grd%im,grd%jm) )         ; grd%scy  = huge(grd%scy(1,1))
  ALLOCATE ( grd%msr(grd%im,grd%jm,grd%km) )  ; grd%msr  = huge(grd%msr(1,1,1))
  ALLOCATE ( grd%imx(grd%km))                 ; grd%imx  = huge(grd%imx(1))
  ALLOCATE (  grd%jmx(grd%km))                ; grd%jmx  = huge(grd%jmx(1))
  ALLOCATE ( grd%istp(grd%im,grd%jm))         ; grd%istp = huge(grd%istp(1,1))
  ALLOCATE ( grd%jstp(grd%im,grd%jm))         ; grd%jstp = huge(grd%jstp(1,1))
  ALLOCATE ( grd%inx(grd%im,grd%jm,grd%km))   ; grd%inx  = huge(grd%inx(1,1,1))
  ALLOCATE ( grd%jnx(grd%im,grd%jm,grd%km))   ; grd%jnx  = huge(grd%jnx(1,1,1))
  ALLOCATE ( grd%fct(grd%im,grd%jm,grd%km) )  ; grd%fct  = huge(grd%fct(1,1,1))
  
  ALLOCATE ( Dump_chl(grd%im,grd%jm,grd%km) ) ; Dump_chl  = 0.0
  ALLOCATE ( Dump_msk(grd%im,grd%jm) )        ; Dump_msk  = 0.0
  ALLOCATE ( grd%chl(grd%im,grd%jm,grd%km,grd%nchl) )    ; grd%chl    = huge(grd%chl(1,1,1,1))
  ALLOCATE ( grd%chl_ad(grd%im,grd%jm,grd%km,grd%nchl) ) ; grd%chl_ad = huge(grd%chl_ad(1,1,1,1))
  
  ALLOCATE ( x3(grd%im,grd%jm,grd%km)) ;  x3 = huge(x3(1,1,1))
  ALLOCATE ( x2(grd%im,grd%jm))        ; x2 = huge(x2(1,1))
  ALLOCATE ( x1(grd%km) )              ;  x1 = huge(x1(1))
  
  if (drv%argo .eq. 1) then
     ALLOCATE ( grd%lon(grd%im,grd%jm)) ; grd%lon = huge(grd%lon(1,1))
     ALLOCATE ( grd%lat(grd%im,grd%jm)) ; grd%lat = huge(grd%lat(1,1))
  endif
  
  ierr = nf90mpi_inq_varid (ncid, 'dx', VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
  ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x2)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all', ierr)
  grd%dx(:,:) = x2(:,:)

  ierr = nf90mpi_inq_varid (ncid, 'dy', VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
  ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x2)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all dy', ierr)
  grd%dy(:,:) = x2(:,:)
  
  if (drv%argo .eq. 1) then
     ierr = nf90mpi_inq_varid (ncid, 'lon', VarId)
     if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
     ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x2)
     if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all lon', ierr)
     grd%lon(:,:) = x2(:,:)

     ierr = nf90mpi_inq_varid (ncid, 'lat', VarId)
     if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
     ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x2)
     if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all lat', ierr)
     grd%lat(:,:) = x2(:,:)
  endif


  ierr = nf90mpi_inq_varid (ncid, 'dep', VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
  ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart(3), MyCount(3), x1)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all dep', ierr)
  grd%dep(:) = x1(:)

  ierr = nf90mpi_inq_varid (ncid, 'tmsk', VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
  ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x3)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all msk', ierr)
  grd%msk(:,:,:) = x3(:,:,:)

  ierr = nf90mpi_inq_varid (ncid, 'regs', VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid regs', ierr)
  ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x2)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all', ierr)
  grd%reg(:,:) = int(x2(:,:))

  ierr = nf90mpi_close(ncid)

  DEALLOCATE ( x3, x2, x1 )


  ! end copy-paste from rdgrds.f90
  ! *****************************************************************************************
  ! *****************************************************************************************

end subroutine parallel_rdgrd

subroutine MyGetDimension(ncid, name, n)
  use pnetcdf
  use mpi
  implicit none

  character name*(*)
  integer :: ncid, ierr
  integer(KIND=MPI_OFFSET_KIND) :: n
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
