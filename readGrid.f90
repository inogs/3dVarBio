subroutine readGrid

  use set_knd
  use drv_str
  use grd_str
  use eof_str
  use filenames
  use mpi_str
  use pnetcdf
  use cns_str
  use bio_str
  use da_params

  implicit none

  integer(i8) :: ierr, ncid, my_km
  integer(i8) :: jpreci, jprecj
  integer(i8) :: VarId
  real(r4), ALLOCATABLE          :: x3(:,:,:), x2(:,:), x1(:)

  integer(8) :: GlobalStart(3), GlobalCount(3)
  integer(KIND=MPI_OFFSET_KIND) MyOffset
  integer    :: MyStatus(MPI_STATUS_SIZE)

  !
  ! open grid1.nc in read-only mode
  ierr = nf90mpi_open(Var3DCommunicator, GRID_FILE, NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_open Grid', ierr)

  !
  ! get grid dimensions
  !
  call MyGetDimension(ncid, 'im', MyOffset)
  GlobalRow = MyOffset

  call MyGetDimension(ncid, 'jm', MyOffset)
  GlobalCol = MyOffset

  call MyGetDimension(ncid, 'km', MyOffset)
  grd%km = MyOffset

  if(grd%km .ne. jpk_200) then
    if(grd%km .gt. jpk_200) then
      if(MyId .eq. 0) then
        write(drv%dia,*) "WARNING!! grd%km differs from jpk_200!!"
        write(drv%dia,*) "grd%km = ", grd%km
        write(drv%dia,*) "jpk_200 = ", jpk_200
        write(drv%dia,*) ""
        write(*,*) "WARNING!! grd%km differs from jpk_200!!"
      endif
    else
      if(MyId .eq. 0) then
        write(drv%dia,*) "Error! grd%km .lt. jpk_200! Aborting"
        write(drv%dia,*) "grd%km = ", grd%km
        write(drv%dia,*) "jpk_200 = ", jpk_200
        write(drv%dia,*) ""
        write(*,*) "Error! grd%km .lt. jpk_200! Aborting"
        write(*,*) "grd%km = ", grd%km
        write(*,*) "jpk_200 = ", jpk_200
        write(*,*) ""
      endif
      call MPI_Barrier(Var3DCommunicator, ierr)
      call MPI_Abort(Var3DCommunicator,-1,ierr)
    endif
  endif

  if(MyId .eq. 0) then
     write(drv%dia,*)'Grid dimensions are: ',GlobalRow,GlobalCol,grd%km

     write(drv%dia,*) ' '
     write(drv%dia,*) 'Dom_size'
     write(drv%dia,*) ' '
     write(drv%dia,*) ' number of processors following i : NumProcI   = ', NumProcI
     write(drv%dia,*) ' number of processors following j : NumProcJ   = ', NumProcJ
     write(drv%dia,*) ' '

     WRITE(*,*) 'Dimension_Med_Grid'
     WRITE(*,*) ' '
     WRITE(*,*) ' GlobalRow  : first  dimension of global domain --> i ',GlobalRow
     WRITE(*,*) ' GlobalCol  : second dimension of global domain --> j ',GlobalCol
     WRITE(*,*) ' Depth      : third dimension of global domain --> j ',grd%km
     WRITE(*,*) ' '
  endif

  GlobalStart(:) = 1
  GlobalCount(1) = GlobalRow
  GlobalCount(2) = GlobalCol
  GlobalCount(3) = grd%km
  ALLOCATE(x3(GlobalRow, GlobalCol, grd%km))
  ALLOCATE(grd%global_msk(GlobalRow, GlobalCol, grd%km))
  ierr = nf90mpi_inq_varid (ncid, 'tmsk', VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
  ierr = nfmpi_get_vara_real_all (ncid, VarId, GlobalStart, GlobalCount, x3)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all global_msk', ierr)
  grd%global_msk(:,:,:) = x3(:,:,:)
  DEALLOCATE(x3)

  call DomainDecomposition


  ALLOCATE(ChlExtended(grd%im+1, grd%jm+1))
  ALLOCATE(SendTop(grd%jm), RecBottom(grd%jm))
  ALLOCATE(SendBottom(grd%jm), RecTop(grd%jm))

  ALLOCATE(SendTop_2d(   grd%jm,grd%km), SendBottom_2d (grd%jm,grd%km))
  ALLOCATE(RecBottom_2d( grd%jm,grd%km), RecTop_2d(     grd%jm,grd%km))
  if(drv%multiv.eq.0) then 
    ALLOCATE(ChlExtended_3d (grd%im+1, grd%jm, grd%km))
  else if(drv%multiv.eq.1) then 
    ALLOCATE(ChlExtended_3d (grd%im+1, grd%jm, ros%kmchl))
  end if
  ALLOCATE(N3nExtended_3d (grd%im+1, grd%jm, grd%km))
  ALLOCATE(O2oExtended_3d (grd%im+1, grd%jm, grd%km))





  ALLOCATE ( grd%reg(grd%im,grd%jm))        ; grd%reg = huge(grd%reg(1,1))
  ALLOCATE ( grd%msk(grd%im,grd%jm,grd%km)) ; grd%msk = huge(grd%msk(1,1,1))
  ALLOCATE ( grd%dep(grd%km))        ; grd%dep = huge(grd%dep(1))
  ALLOCATE ( grd%dx(grd%im,grd%jm))  ; grd%dx  = huge(grd%dx(1,1))
  ALLOCATE ( grd%dy(grd%im,grd%jm))  ; grd%dy  = huge(grd%dy(1,1))

  ALLOCATE ( grd%alx(GlobalRow,localCol,grd%km) )         ; grd%alx  = huge(grd%alx(1,1,1))
  ALLOCATE ( grd%aly(localRow,GlobalCol,grd%km) )         ; grd%aly  = huge(grd%aly(1,1,1))
  ALLOCATE ( grd%btx(GlobalRow,localCol,grd%km) )         ; grd%btx  = huge(grd%btx(1,1,1))
  ALLOCATE ( grd%bty(localRow,GlobalCol,grd%km) )         ; grd%bty  = huge(grd%bty(1,1,1))
  ALLOCATE ( grd%scx(GlobalRow,localCol,grd%km) )         ; grd%scx  = huge(grd%scx(1,1,1))
  ALLOCATE ( grd%scy(localRow,GlobalCol,grd%km) )         ; grd%scy  = huge(grd%scy(1,1,1))
  ALLOCATE ( grd%imx(grd%km))                 ; grd%imx  = huge(grd%imx(1))
  ALLOCATE (  grd%jmx(grd%km))                ; grd%jmx  = huge(grd%jmx(1))
  ALLOCATE ( grd%istp(GlobalRow,localCol,grd%km))         ; grd%istp = huge(grd%istp(1,1,1))
  ALLOCATE ( grd%jstp(localRow,GlobalCol,grd%km))         ; grd%jstp = huge(grd%jstp(1,1,1))
  ALLOCATE ( grd%inx(GlobalRow,localCol,grd%km))   ; grd%inx  = huge(grd%inx(1,1,1))
  ALLOCATE ( grd%jnx(localRow,GlobalCol,grd%km))   ; grd%jnx  = huge(grd%jnx(1,1,1))

  if(drv%multiv .eq. 0) then
    if(drv%chl_assim .eq. 1) then
      ALLOCATE ( grd%chl(grd%im,grd%jm,grd%km) )    ; grd%chl    = huge(grd%chl(1,1,1))
      ALLOCATE ( grd%chl_ad(grd%im,grd%jm,grd%km) ) ; grd%chl_ad = huge(grd%chl_ad(1,1,1))
      ALLOCATE ( bio%phy(grd%im,grd%jm,grd%km,bio%nphy,bio%ncmp) ) ; bio%phy = huge(bio%phy(1,1,1,1,1))
      ALLOCATE ( bio%phy_ad(grd%im,grd%jm,grd%km,bio%nphy,bio%ncmp) ) ; bio%phy_ad = huge(bio%phy_ad(1,1,1,1,1))
    endif
    if(drv%nut .eq. 1) then
      if(bio%N3n .eq. 1) then
        ALLOCATE ( grd%n3n(grd%im,grd%jm,grd%km) )    ; grd%n3n    = huge(grd%n3n(1,1,1))
        ALLOCATE ( grd%n3n_ad(grd%im,grd%jm,grd%km) ) ; grd%n3n_ad = huge(grd%n3n_ad(1,1,1))
      endif
      if(bio%O2o .eq. 1) then
        ALLOCATE ( grd%o2o(grd%im,grd%jm,grd%km) )    ; grd%o2o    = huge(grd%o2o(1,1,1))
        ALLOCATE ( grd%o2o_ad(grd%im,grd%jm,grd%km) ) ; grd%o2o_ad = huge(grd%o2o_ad(1,1,1))
      endif
    endif

  else if(drv%multiv .eq. 1) then
    ALLOCATE ( grd%chl(grd%im,grd%jm,ros%kmchl) )    ; grd%chl    = huge(grd%chl(1,1,1))
    ALLOCATE ( grd%chl_ad(grd%im,grd%jm,ros%kmchl) ) ; grd%chl_ad = huge(grd%chl_ad(1,1,1))
    ALLOCATE ( grd%n3n(grd%im,grd%jm,ros%kmnit) )    ; grd%n3n    = huge(grd%n3n(1,1,1))
    ALLOCATE ( grd%n3n_ad(grd%im,grd%jm,ros%kmnit) ) ; grd%n3n_ad = huge(grd%n3n_ad(1,1,1))
    ALLOCATE ( bio%phy(grd%im,grd%jm,ros%kmchl,bio%nphy,bio%ncmp) ) ; bio%phy = huge(bio%phy(1,1,1,1,1))
    ALLOCATE ( bio%phy_ad(grd%im,grd%jm,ros%kmchl,bio%nphy,bio%ncmp) ) ; bio%phy_ad = huge(bio%phy_ad(1,1,1,1,1))
  endif

  
  ALLOCATE ( x3(grd%im,grd%jm,grd%km)) ;  x3 = huge(x3(1,1,1))
  ALLOCATE ( x2(grd%im,grd%jm))        ;  x2 = huge(x2(1,1))
  ALLOCATE ( x1(grd%km) )              ;  x1 = huge(x1(1))

  if (drv%argo_obs .eq. 1) then
     ALLOCATE ( grd%lon(grd%im,grd%jm)) ; grd%lon = huge(grd%lon(1,1))
     ALLOCATE ( grd%lat(grd%im,grd%jm)) ; grd%lat = huge(grd%lat(1,1))
  endif

  ierr = nf90mpi_inq_varid (ncid, 'dx', VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
  ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x2)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all dx', ierr)
  grd%dx(:,:) = x2(:,:)

  ierr = nf90mpi_inq_varid (ncid, 'dy', VarId)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
  ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x2)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all dy', ierr)
  grd%dy(:,:) = x2(:,:)

  if (drv%argo_obs .eq. 1) then
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
     
     grd%NextLongitude=grd%lon(1,1)
     ! Send to ProcTop with Tag = MyId and receiving from 
     ! ProcBottom with Tag = ProcBottom :)
     call MPI_Sendrecv_replace(grd%NextLongitude,1,MPI_REAL8,ProcTop,MyId,&
      ProcBottom,ProcBottom, Var3DCommunicator, MyStatus, ierr)
     if(ProcBottom .eq. MPI_PROC_NULL) grd%NextLongitude = grd%lon(grd%im,grd%jm)
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
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all regs', ierr)
  grd%reg(:,:) = int(x2(:,:))

  ierr = nf90mpi_close(ncid)

  DEALLOCATE ( x3, x2, x1 )


  ! end copy-paste from rdgrds.f90
  ! *****************************************************************************************
  ! *****************************************************************************************


end subroutine readGrid

subroutine DomainDecomposition

  use drv_str
  use mpi_str
  use grd_str
  use eof_str

  implicit none

  integer, allocatable :: ilcit(:,:), ilcjt(:,:), BalancedSlice(:,:)
  integer(i8) :: ji, jj, TmpInt, ierr ! jpi, jpj, nn, i
  integer(i8) :: GlobalRestCol, GlobalRestRow
  integer(i8) :: i, j, k, kk, my_km
  integer(i8) :: NCoastX, NCoastY, TmpCoast
  integer(i8) :: NRows, NCols
  integer(i8) :: SliceRestRow, SliceRestCol
  integer(i8) :: OffsetCol, OffsetRow
  real        :: TotX, TotY, C
  integer     :: nnx, ii, iProc, i0

  integer, allocatable :: ToBalanceX(:,:), ToBalanceY(:,:)
  integer, allocatable, dimension(:) :: SendDisplY2D, SendCountY2D

  GlobalRestRow = mod(GlobalRow, NumProcI)
  GlobalRestCol = mod(GlobalCol, NumProcJ)

  !*******************************************
  !
  ! PDICERBO version of the domain decomposition:
  ! the domain is divided among the processes into slices
  ! of size (GlobalRow / NumProcI, GlobalCol / NumProcJ).
  ! Clearly, the division is done tacking into account
  ! rests. The only condition we need is that NumProcI*NumProcJ = NPROC
  !
  ! WARNING!!! netcdf stores data in ROW MAJOR order
  ! while here we are reading in column major order.
  ! We have to take into account this simply swapping ("ideally")
  ! the entries of MyStart and MyCount
  !
  !*******************************************

  ! computing rests for X direction
  MyCount(1) = GlobalRow / NumProcI
  MyCount(2) = GlobalCol / NumProcJ

  OffsetRow = 0
  if (MyPosI .lt. GlobalRestRow) then
    MyCount(1) = MyCount(1) + 1
    OffsetRow = MyPosI
  else
    OffsetRow = GlobalRestRow
  end if

  if(MyPosI+1 .eq. GlobalRestRow) then
    NextLocalRow = MyCount(1)-1
  else
    NextLocalRow = MyCount(1)
  endif

  ! computing rests for Y direction
  OffsetCol = 0
  if (MyPosJ .lt. GlobalRestCol) then
    MyCount(2) = MyCount(2) + 1
  else
    OffsetCol = GlobalRestCol
  end if

  TmpInt = GlobalRow / NumProcI
  MyStart(1) = TmpInt * MyPosI + OffsetRow + 1
  MyCount(1) = MyCount(1)

  TmpInt = MyId / NumProcI
  MyStart(2) = mod(MyCount(2) * TmpInt + OffsetCol, GlobalCol) + 1
  MyCount(2) = MyCount(2)

  ! taking all values along k direction
  MyStart(3) = 1
  MyCount(3) = grd%km

  if(drv%Verbose .eq. 1) &
       write(*,*) "MyId = ", MyId, " MyStart = ", MyStart, " MyCount = ", MyCount

  grd%im = MyCount(1)
  grd%jm = MyCount(2)

  !
  ! initializing quantities needed to slicing along i and j directions
  !
  localRow = grd%im / NumProcJ
  localCol = grd%jm / NumProcI
  SliceRestRow = mod(grd%im, NumProcJ)
  SliceRestCol = mod(grd%jm, NumProcI)

  ! x direction (-> GlobalRow)
  if(SliceRestCol .ne. 0) then
    if(MyPosI .lt. SliceRestCol) &
          localCol = localCol + 1
  end if

  SendDisplX3D(1) = 0
  RecDisplX3D(1)  = 0
  SendDisplX3D_chl(1) = 0
  RecDisplX3D_chl(1)  = 0

  SendDisplX2D(1) = 0
  RecDisplX2D(1)  = 0

  my_km = grd%km
  if(drv%multiv.eq.1) &
    my_km = ros%kmchl

  do i=1,NumProcI
    if(i-1 .lt. SliceRestCol) then
        OffsetRow = 1
    else
        OffsetRow = 0
    end if

    if(i-1 .lt. mod(GlobalRow, NumProcI)) then
        OffsetCol = 1
    else
        OffsetCol = 0
    end if

    SendCountX3D(i) = (grd%jm / NumProcI + OffsetRow) * grd%im * grd%km
    RecCountX3D(i)  = localCol * grd%km * (GlobalRow / NumProcI + OffsetCol)

    SendCountX3D_chl(i) = (grd%jm / NumProcI + OffsetRow) * grd%im * my_km
    RecCountX3D_chl(i)  = localCol * my_km * (GlobalRow / NumProcI + OffsetCol)

    SendCountX2D(i) = (grd%jm / NumProcI + OffsetRow) * grd%im
    RecCountX2D(i)  = localCol * (GlobalRow / NumProcI + OffsetCol)

    if(i .lt. NumProcI) then
        SendDisplX3D(i+1) = SendDisplX3D(i) + SendCountX3D(i)
        RecDisplX3D(i+1)  = RecDisplX3D(i) + RecCountX3D(i)

        SendDisplX3D_chl(i+1) = SendDisplX3D_chl(i) + SendCountX3D_chl(i)
        RecDisplX3D_chl(i+1)  = RecDisplX3D_chl(i) + RecCountX3D_chl(i)

        SendDisplX2D(i+1) = SendDisplX2D(i) + SendCountX2D(i)
        RecDisplX2D(i+1)  = RecDisplX2D(i) + RecCountX2D(i)
    end if
  end do

  ALLOCATE(SendDisplY2D(NumProcJ), SendCountY2D(NumProcJ))
  ! y direction (-> GlobalCol)
  if(SliceRestRow .ne. 0) then
    if(MyPosJ .lt. SliceRestRow) &
          localRow = localRow + 1
  end if

  SendDisplY2D(1) = 0

  do i=1,NumProcJ
    if(i-1 .lt. SliceRestRow) then
        OffsetCol = 1
    else
        OffsetCol = 0
    end if

    if(i-1 .lt. mod(GlobalCol, NumProcJ)) then
        OffsetRow = 1
    else
        OffsetRow = 0
    end if

    SendCountY2D(i) = (grd%im / NumProcJ + OffsetCol) * grd%jm

    if(i .lt. NumProcJ) then
        SendDisplY2D(i+1) = SendDisplY2D(i) + SendCountY2D(i)
    end if
  end do

  if(MyPosI .lt. GlobalRestRow) then
    TmpInt = 0
  else
    TmpInt = 1
  end if
  GlobalRowOffset = SendDisplY2D(MyPosJ+1)/grd%jm + MyPosI*grd%im + TmpInt*GlobalRestRow

  if(MyPosJ .lt. GlobalRestCol) then
    TmpInt = 0
  else
    TmpInt = 1
  end if
  GlobalColOffset = SendDisplX2D(MyPosI+1)/grd%im + MyPosJ*grd%jm + TmpInt*GlobalRestCol

  DEALLOCATE(SendDisplY2D, SendCountY2D)

end subroutine DomainDecomposition

subroutine MyMax(arr, GlobCol, km, i0, ii, k, val)

  implicit none

  integer :: i0, ii, k, j, GlobCol, km
  integer :: arr(GlobCol,km)
  integer, intent(inout)    :: val

  do j=i0,ii
     val = max(val, arr(j, k))
  end do

end subroutine MyMax

subroutine MyGetDimension(ncid, name, n)
  use pnetcdf
  use mpi_str
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

  use mpi_str
  use pnetcdf

  implicit none

  character*(*), intent(in) :: err_msg
  integer,       intent(in) :: errcode

  !local variables
  integer err

  write(*,*) 'Error: ', trim(err_msg), ' ', nf90mpi_strerror(errcode)
  call MPI_Abort(Var3DCommunicator, -1, err)
  return
end subroutine handle_err



