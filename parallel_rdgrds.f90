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
  integer(i8) :: TmpInt, VarId
  real(r4), ALLOCATABLE          :: x3(:,:,:), x2(:,:), x1(:)

  integer, allocatable :: ilcit(:,:), ilcjt(:,:)
  integer(i8) :: ji, jj ! jpi, jpj, nn, i 
  integer(i8) :: GlobalRestCol, GlobalRestRow, i
  integer(i8) :: SliceRestRow, SliceRestCol
  integer(i8) :: OffsetCol, OffsetRow

  integer(8) :: GlobalStart(3), GlobalCount(3)
  integer(KIND=MPI_OFFSET_KIND) MyOffset
  
  !
  ! open grid1.nc in read-only mode
  ierr = nf90mpi_open(MPI_COMM_WORLD, GRID_FILE, NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_open', ierr)

  !
  ! get grid dimensions
  !
  call MyGetDimension(ncid, 'im', MyOffset)
  GlobalRow = MyOffset

  call MyGetDimension(ncid, 'jm', MyOffset)
  GlobalCol = MyOffset

  call MyGetDimension(ncid, 'km', MyOffset)
  grd%km = MyOffset

  if(MyRank .eq. 0) then
     write(drv%dia,*)'Grid dimensions are: ',GlobalRow,GlobalCol,grd%km

     WRITE(*,*) 'Dimension_Med_Grid'
     WRITE(*,*) ' '
     WRITE(*,*) ' GlobalRow  : first  dimension of global domain --> i ',GlobalRow
     WRITE(*,*) ' GlobalCol  : second dimension of global domain --> j ',GlobalCol
     WRITE(*,*) ' ReadDomDec : ',drv%ReadDomDec
     WRITE(*,*) ' '
  endif
  
  GlobalRestRow = mod(GlobalRow, NumProcI)
  GlobalRestCol = mod(GlobalCol, NumProcJ)
  
  if(drv%ReadDomDec .eq. 1) then
     allocate(ilcit(NumProcI, NumProcJ)) ; ilcit = huge(ilcit(1,1))
     allocate(ilcjt(NumProcI, NumProcJ)) ; ilcjt = huge(ilcjt(1,1))
     
     open(3333,file='Dom_Dec_jpi.ascii', form='formatted')
     open(3334,file='Dom_Dec_jpj.ascii', form='formatted')
     
     read(3333,*) ((ilcit(ji,jj), jj=1,NumProcJ),ji=1,NumProcI)
     read(3334,*) ((ilcjt(ji,jj), jj=1,NumProcJ),ji=1,NumProcI)
  
     close(3333)
     close(3334)

     do ji=1,NumProcI
        do jj=1,NumProcJ
           if(NumProcJ .gt. 1) then
              if(mod(NumProcJ-jj, NumProcJ-1) .eq. 0) then
                 ilcjt(ji,jj) = ilcjt(ji,jj) - 1
              else
                 ilcjt(ji,jj) = ilcjt(ji,jj) - 2
              end if
           end if
           if(NumProcI .gt. 1) then
              if(mod(NumProcI-ji, NumProcI-1) .eq. 0) then
                 ilcit(ji,jj) = ilcit(ji,jj) - 1
              else
                 ilcit(ji,jj) = ilcit(ji,jj) - 2
              end if
           end if
        end do
     end do

     grd%im = ilcit(MyPosI+1, MyPosJ+1)
     grd%jm = ilcjt(MyPosI+1, MyPosJ+1)
     MyCount(1) = grd%im
     MyCount(2) = grd%jm
     MyCount(3) = grd%km

     MyStart(:) = 1

     do i=1,MyPosI
        MyStart(1) = MyStart(1) + ilcit(i,MyPosJ+1)
     end do
     do i=1,MyPosJ
        MyStart(2) = MyStart(2) + ilcjt(MyPosI+1, i)
     end do
     MyStart(3) = 1

     write(*,*) "MyRank = ", MyRank, " MyStart = ", MyStart, " MyCount = ", MyCount

     !
     ! initializing quantities needed to slicing along i and j directions
     !
     localRow = grd%im / NumProcJ
     localCol = grd%jm / NumProcI
     SliceRestRow = mod(grd%im, NumProcJ)
     SliceRestCol = mod(grd%jm, NumProcI)

     ! x direction (-> GlobalRow)
     if(SliceRestCol .ne. 0) then !print*,"WARNING!!!!!! mod(grd%jm, NumProcI) .ne. 0!!! Case not implemented yet!!"
        if(MyPosI .lt. SliceRestCol) &
             localCol = localCol + 1
     end if
     
     SendDisplX4D(1) = 0
     RecDisplX4D(1)  = 0
     
     SendDisplX2D(1) = 0
     RecDisplX2D(1)  = 0
     
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
        
        SendCountX4D(i) = (grd%jm / NumProcI + OffsetRow) * grd%im * grd%km
        RecCountX4D(i)  = localCol * grd%km * ilcit(i, MyPosJ+1) ! (GlobalRow / NumProcI + OffsetCol)
        
        SendCountX2D(i) = (grd%jm / NumProcI + OffsetRow) * grd%im
        RecCountX2D(i)  = localCol * ilcit(i, MyPosJ+1) ! (GlobalRow / NumProcI + OffsetCol)
        
        if(i .lt. NumProcI) then
           SendDisplX4D(i+1) = SendDisplX4D(i) + SendCountX4D(i)
           RecDisplX4D(i+1)  = RecDisplX4D(i) + RecCountX4D(i)
           
           SendDisplX2D(i+1) = SendDisplX2D(i) + SendCountX2D(i)
           RecDisplX2D(i+1)  = RecDisplX2D(i) + RecCountX2D(i)
        end if
     end do

     ! y direction (-> GlobalCol)
     if(SliceRestRow .ne. 0) then
        if(MyPosJ .lt. SliceRestRow) &
             localRow = localRow + 1
     end if
     
     SendDisplY4D(1) = 0
     RecDisplY4D(1)  = 0
     
     SendDisplY2D(1) = 0
     RecDisplY2D(1)  = 0
     
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
        
        SendCountY4D(i) = (grd%im / NumProcJ + OffsetCol) * grd%jm * grd%km
        RecCountY4D(i)  = localRow * grd%km * ilcjt(MyPosI+1, i)
        
        SendCountY2D(i) = (grd%im / NumProcJ + OffsetCol) * grd%jm
        RecCountY2D(i)  = localRow * ilcjt(MyPosI+1, i)
        
        if(i .lt. NumProcJ) then
           SendDisplY4D(i+1) = SendDisplY4D(i) + SendCountY4D(i)
           RecDisplY4D(i+1)  = RecDisplY4D(i) + RecCountY4D(i)
           
           SendDisplY2D(i+1) = SendDisplY2D(i) + SendCountY2D(i)
           RecDisplY2D(i+1)  = RecDisplY2D(i) + RecCountY2D(i)
        end if
     end do

     if(MyPosI .lt. GlobalRestRow) then
        TmpInt = 0
     else
        TmpInt = 1
     end if
     GlobalRowOffset = SendDisplY2D(MyPosJ+1)/grd%jm
     do i=1,MyPosI
        GlobalRowOffset = GlobalRowOffset + ilcit(i, MyPosJ+1)
     end do
  
     if(MyPosJ .lt. GlobalRestCol) then
        TmpInt = 0
     else
        TmpInt = 1
     end if
     GlobalColOffset = SendDisplX2D(MyPosI+1)/grd%im
     do i=1,MyPosJ
        GlobalColOffset = GlobalColOffset + ilcjt(MyPosI+1, i)
     end do
  
  else ! drv%ReadDomDec .eq. 0
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
     
     TmpInt = MyRank / NumProcI
     MyStart(2) = mod(MyCount(2) * TmpInt + OffsetCol, GlobalCol) + 1
     MyCount(2) = MyCount(2)
     
     ! taking all values along k direction
     MyStart(3) = 1
     MyCount(3) = grd%km
     
     write(*,*) "MyRank = ", MyRank, " MyStart = ", MyStart, " MyCount = ", MyCount
     
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
     
     SendDisplX4D(1) = 0
     RecDisplX4D(1)  = 0
     
     SendDisplX2D(1) = 0
     RecDisplX2D(1)  = 0
     
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
        
        SendCountX4D(i) = (grd%jm / NumProcI + OffsetRow) * grd%im * grd%km
        RecCountX4D(i)  = localCol * grd%km * (GlobalRow / NumProcI + OffsetCol)
        
        SendCountX2D(i) = (grd%jm / NumProcI + OffsetRow) * grd%im
        RecCountX2D(i)  = localCol * (GlobalRow / NumProcI + OffsetCol)
        
        if(i .lt. NumProcI) then
           SendDisplX4D(i+1) = SendDisplX4D(i) + SendCountX4D(i)
           RecDisplX4D(i+1)  = RecDisplX4D(i) + RecCountX4D(i)
           
           SendDisplX2D(i+1) = SendDisplX2D(i) + SendCountX2D(i)
           RecDisplX2D(i+1)  = RecDisplX2D(i) + RecCountX2D(i)
        end if
     end do
     
     ! y direction (-> GlobalCol)
     if(SliceRestRow .ne. 0) then
        if(MyPosJ .lt. SliceRestRow) &
             localRow = localRow + 1
     end if
     
     SendDisplY4D(1) = 0
     RecDisplY4D(1)  = 0
     
     SendDisplY2D(1) = 0
     RecDisplY2D(1)  = 0
     
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
        
        SendCountY4D(i) = (grd%im / NumProcJ + OffsetCol) * grd%jm * grd%km
        RecCountY4D(i)  = localRow * grd%km * (GlobalCol / NumProcJ + OffsetRow)
        
        SendCountY2D(i) = (grd%im / NumProcJ + OffsetCol) * grd%jm
        RecCountY2D(i)  = localRow * (GlobalCol / NumProcJ + OffsetRow)
        
        if(i .lt. NumProcJ) then
           SendDisplY4D(i+1) = SendDisplY4D(i) + SendCountY4D(i)
           RecDisplY4D(i+1)  = RecDisplY4D(i) + RecCountY4D(i)
           
           SendDisplY2D(i+1) = SendDisplY2D(i) + SendCountY2D(i)
           RecDisplY2D(i+1)  = RecDisplY2D(i) + RecCountY2D(i)
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
     
  end if ! drv%ReadDomDec .eq. 0
  
  ALLOCATE(ChlExtended(grd%im+1, grd%jm+1, grd%nchl))
  ALLOCATE(SendLeft(grd%im), RecRight(grd%im))
  ALLOCATE(SendRight(grd%im), RecLeft(grd%im))
  ALLOCATE(SendTop(grd%jm), RecBottom(grd%jm))
  ALLOCATE(SendBottom(grd%jm), RecTop(grd%jm))  
  
  ALLOCATE(ChlExtendedAD_4D(0:(grd%im+1), 0:(grd%jm+1), grd%km, grd%nchl))
  ALLOCATE(ChlExtended4D(0:(grd%im+1), 0:(grd%jm+1), grd%km, grd%nchl))
  ALLOCATE(SendLeft2D(grd%im, grd%km), RecRight2D(grd%im, grd%km))
  ALLOCATE(SendRight2D(grd%im, grd%km), RecLeft2D(grd%im, grd%km))
  ALLOCATE(SendTop2D(grd%jm, grd%km), RecBottom2D(grd%jm, grd%km))
  ALLOCATE(SendBottom2D(grd%jm, grd%km), RecTop2D(grd%jm, grd%km))
  
  ! print*, "Debugging", MyRank, "SC", SendCountY2D, "RC", RecCountY2D, "SD", SendDisplY2D, "RD", RecDisplY2D
  
  ! *****************************************************************************************
  ! *****************************************************************************************
  ! (almost) copy-paste from rdgrds.f90
  ! Allocate grid arrays

  ALLOCATE ( grd%reg(grd%im,grd%jm))        ; grd%reg = huge(grd%reg(1,1))
  ALLOCATE ( grd%msk(grd%im,grd%jm,grd%km)) ; grd%msk = huge(grd%msk(1,1,1))
  ALLOCATE ( grd%dep(grd%km))        ; grd%dep = huge(grd%dep(1))
  ALLOCATE ( grd%dx(grd%im,grd%jm))  ; grd%dx  = huge(grd%dx(1,1))
  ALLOCATE ( grd%dy(grd%im,grd%jm))  ; grd%dy  = huge(grd%dy(1,1))
  
  ALLOCATE ( grd%alx(GlobalRow,localCol) )         ; grd%alx  = huge(grd%alx(1,1))
  ALLOCATE ( grd%aly(localRow,GlobalCol) )         ; grd%aly  = huge(grd%aly(1,1))
  ALLOCATE ( grd%btx(GlobalRow,localCol) )         ; grd%btx  = huge(grd%btx(1,1))
  ALLOCATE ( grd%bty(localRow,GlobalCol) )         ; grd%bty  = huge(grd%bty(1,1))
  ALLOCATE ( grd%scx(GlobalRow,localCol) )         ; grd%scx  = huge(grd%scx(1,1))
  ALLOCATE ( grd%scy(localRow,GlobalCol) )         ; grd%scy  = huge(grd%scy(1,1))
  ALLOCATE ( grd%msr(grd%im,grd%jm,grd%km) )  ; grd%msr  = huge(grd%msr(1,1,1))
  ALLOCATE ( grd%imx(grd%km))                 ; grd%imx  = huge(grd%imx(1))
  ALLOCATE (  grd%jmx(grd%km))                ; grd%jmx  = huge(grd%jmx(1))
  ALLOCATE ( grd%istp(GlobalRow,localCol))         ; grd%istp = huge(grd%istp(1,1))
  ALLOCATE ( grd%jstp(localRow,GlobalCol))         ; grd%jstp = huge(grd%jstp(1,1))
  ALLOCATE ( grd%inx(GlobalRow,localCol,grd%km))   ; grd%inx  = huge(grd%inx(1,1,1))
  ALLOCATE ( grd%jnx(localRow,GlobalCol,grd%km))   ; grd%jnx  = huge(grd%jnx(1,1,1))
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

  GlobalStart(:) = 1
  GlobalCount(1) = GlobalRow
  GlobalCount(2) = GlobalCol
  GlobalCount(3) = grd%km
  DEALLOCATE(x3)
  ALLOCATE(grd%global_msk(GlobalRow, GlobalCol, grd%km))
  ALLOCATE(x3(GlobalRow, GlobalCol, grd%km))
  ierr = nfmpi_get_vara_real_all (ncid, VarId, GlobalStart, GlobalCount, x3)
  if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all global_msk', ierr)
  grd%global_msk(:,:,:) = x3(:,:,:)
  
  
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
