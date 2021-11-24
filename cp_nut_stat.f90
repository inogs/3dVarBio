subroutine cp_nut_stat

  use set_knd
  use grd_str
  use drv_str
  use mpi_str
  use bio_str
  use pnetcdf
  use da_params

  implicit none

  INTEGER(i4)        :: ncid, ierr, i, j, k, l
  INTEGER(i4)        :: idP, iVar
  INTEGER(I4)        :: xid,yid,depid,timeId, idTim
  INTEGER            :: system, SysErr

  INTEGER(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km, MyTime
  INTEGER(KIND=MPI_OFFSET_KIND) :: MyCountSingle(1), MyStartSingle(1)
  CHARACTER(LEN=46)    :: BioRestart
  CHARACTER(LEN=47)    :: BioRestartLong
  CHARACTER(LEN=6)     :: MyVarName
  ! LOGICAL, ALLOCATABLE :: MyConditions(:,:,:,:)

  ! bug fix Intel 2018
  real(r4), allocatable, dimension(:,:,:,:) :: DumpBio
  integer(KIND=MPI_OFFSET_KIND) :: MyStart_4d(4), MyCount_4d(4)

  real(r8) :: TimeArr(1)


  MyStart_4d(1:3) = MyStart(:)
  MyStart_4d(4) = 1
  MyCount_4d(1:3) = MyCount(:)
  MyCount_4d(4) = 1

  
  ALLOCATE(DumpBio(grd%im,grd%jm,grd%km,1)); DumpBio(:,:,:,:) = 1.e20
  ! ALLOCATE(MyConditions(grd%im,grd%jm,grd%km,bio%nphy))

  if(MyId .eq. 0) then
     write(drv%dia,*) 'writing nut structure (only copy from RSTbefore)'     
     write(*,*) 'writing nut structure (only copy from RSTbefore)'          
  endif

  global_im = GlobalRow
  global_jm = GlobalCol
  global_km = grd%km
  MyTime = 1

  MyCountSingle(1) = 1
  MyStartSingle(1) = 1
  TimeArr(1) = DA_JulianDate


  do l=1,NNutVar
    iVar = NPhytoVar + l

    if(iVar .gt. NBioVar) then
      if(MyId .eq. 0) &
        write(*,*) "Warning: Reading a variable not in the DA_VarList!"
    endif

    BioRestart = 'DA__FREQ_1/RST_after.'//ShortDate//'.'//DA_VarList(iVar)//'.nc'
    BioRestartLong = 'DA__FREQ_1/RST_after.'//DA_DATE//'.'//DA_VarList(iVar)//'.nc'

    if(drv%Verbose .eq. 1 .and. MyId .eq. 0) &
      print*, "Writing Nut Restart ", BioRestart
    
    ierr = nf90mpi_create(Var3DCommunicator, BioRestart, NF90_CLOBBER, MPI_INFO_NULL, ncid)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_create '//BioRestart, ierr)

    ierr = nf90mpi_def_dim(ncid,'x',global_im ,xid)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim longitude ', ierr)
    ierr = nf90mpi_def_dim(ncid,'y' ,global_jm ,yid)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim latitude ', ierr)
    ierr = nf90mpi_def_dim(ncid,'z'    ,global_km, depid)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim depth ', ierr)
    ierr = nf90mpi_def_dim(ncid,'time',MyTime ,timeId)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim time ', ierr)

    MyVarName='TRN'//DA_VarList(iVar)

    ierr = nf90mpi_def_var(ncid, MyVarName, nf90_float, (/xid,yid,depid,timeId/), idP )
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)
    
    ierr = nf90mpi_def_var(ncid,'time'   , nf90_double,  (/timeid/)  , idTim)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)
    ierr = nf90mpi_put_att(ncid,idP   , 'missing_value',1.e+20)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_att', ierr)

    ierr = nf90mpi_enddef(ncid)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_enddef'//DA_VarList(iVar), ierr)

    do k=1,grd%km
      do j=1,grd%jm
        do i=1,grd%im

          if(bio%InitialNut(i,j,k,1) .lt. 1.e20) then
              DumpBio(i,j,k,1) = bio%InitialNut(i,j,k,l)
          endif

        enddo
      enddo
    enddo

    ierr = nf90mpi_put_var_all(ncid,idP,DumpBio,MyStart_4d,MyCount_4d)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//DA_VarList(iVar), ierr)

    ierr = nf90mpi_put_var_all(ncid,idTim,TimeArr,MyStartSingle,MyCountSingle)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//DA_VarList(iVar), ierr)

    ierr = nf90mpi_close(ncid)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_close '//BioRestart, ierr)

    call MPI_Barrier(Var3DCommunicator, ierr)
    ! only process 0 creates link to restart files
    if(MyId .eq. 0) then
      SysErr = system("ln -sf $PWD/"//BioRestart//" "//BioRestartLong)
      if(SysErr /= 0) call MPI_Abort(MPI_COMM_WORLD, -1, SysErr)
    endif
  enddo ! l

  DEALLOCATE(DumpBio)
  ! DEALLOCATE(DumpBio, ValuesToTest, MyConditions)

end subroutine cp_nut_stat
