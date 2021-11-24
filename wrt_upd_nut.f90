subroutine wrt_upd_nut

  use set_knd
  use grd_str
  use drv_str
  use mpi_str
  use bio_str
  use pnetcdf
  use da_params

  implicit none

  INTEGER(i4)        :: ncid, ierr, i, j, k, l, m
  INTEGER(i4)        :: idP, iVar
  INTEGER(I4)        :: xid,yid,depid,timeId, idTim
  INTEGER            :: system, SysErr

  INTEGER(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km, MyTime
  INTEGER(KIND=MPI_OFFSET_KIND) :: MyCountSingle(1), MyStartSingle(1)
  CHARACTER(LEN=45)    :: BioRestart
  CHARACTER(LEN=47)    :: BioRestartLong
  CHARACTER(LEN=6)     :: MyVarName


  real(r8)           :: TmpVal, Myn3nUpdate, Myn1pUpdate, totchl, SMALL
  real(r4), allocatable, dimension(:,:,:,:) :: ValuesToTest

  ! bug fix Intel 2018
  real(r4), allocatable, dimension(:,:,:,:) :: DumpBio
  integer(KIND=MPI_OFFSET_KIND) :: MyStart_4d(4), MyCount_4d(4)

  real(r8) :: TimeArr(1)

  SMALL     =  1.e-5

  MyStart_4d(1:3) = MyStart(:)
  MyStart_4d(4) = 1
  MyCount_4d(1:3) = MyCount(:)
  MyCount_4d(4) = 1

  
  ALLOCATE(DumpBio(grd%im,grd%jm,grd%km,1)); DumpBio(:,:,:,:) = 1.e20
  ALLOCATE(ValuesToTest(grd%im,grd%jm,grd%km,2)); ValuesToTest(:,:,:,:) = dble(0.)

  if(MyId .eq. 0) then
     write(drv%dia,*) 'writing  nut update based on chl'
     write(*,*) 'writing  nut update based on chl'
  endif

  global_im = GlobalRow
  global_jm = GlobalCol
  global_km = grd%km
  MyTime = 1

  MyCountSingle(1) = 1
  MyStartSingle(1) = 1
  TimeArr(1) = DA_JulianDate


  do k=1,grd%km
    do j=1,grd%jm
      do i=1,grd%im
        totchl = 0.
        do l=1,bio%nphy
           totchl = totchl + bio%phy(i,j,k,l,1)
        enddo
        Myn3nUpdate = totchl*bio%covn3n_chl(i,j,k)
        ValuesToTest(i,j,k,1) = bio%InitialNut(i,j,k,1) + Myn3nUpdate

        Myn1pUpdate = totchl*bio%covn1p_chl(i,j,k)
        ValuesToTest(i,j,k,2) = bio%InitialNut(i,j,k,2) + Myn1pUpdate

      enddo
    enddo
  enddo

  ! do k=1,grd%km
  !   do j=1,grd%jm
  !     do i=1,grd%im
  !       totchl = 0.
  !       do l=1,bio%nphy
  !          totchl = totchl + bio%phy(i,j,k,l,1)
  !       enddo
  !       Myn3nUpdate = totchl*bio%covn3n_chl(i,j,k)
  !       ValuesToTest(i,j,k,1) = bio%InitialNut(i,j,k,1) + Myn3nUpdate

  !       ValuesToTest(i,j,k,2) = bio%InitialNut(i,j,k,2) + Myn3nUpdate*bio%covn3n_n1p(i,j,k)

  !     enddo
  !   enddo
  ! enddo

  do l=1,NNutVar
    iVar = NPhytoVar + l

    if(iVar .gt. NBioVar) then
      if(MyId .eq. 0) &
        write(*,*) "Warning: Reading a variable not in the DA_VarList!"
    endif


    BioRestart = 'DA__FREQ_1/RST_after.'//ShortDate//'.'//DA_VarList(iVar)//'.nc'
    BioRestartLong = 'DA__FREQ_1/RST_after.'//DA_DATE//'.'//DA_VarList(iVar)//'.nc'

    if(drv%Verbose .eq. 1 .and. MyId .eq. 0) &
      print*, "Writing Nut Restart based on chl ", BioRestart
      
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

          if(bio%InitialNut(i,j,k,l) .lt. 1.e20) then

            if(grd%msk(i,j,k).eq.1) then

              if(ValuesToTest(i,j,k,l) .lt. 0) then
                  ! Excluding negative concentrations
                  ! This correction must be the first
                  ! condition applied (before apply corrections
                  ! on the other components)
                TmpVal = 0.01*bio%InitialNut(i,j,k,l)
                if(TmpVal.gt.SMALL) then
                  TmpVal = SMALL
                endif
                DumpBio(i,j,k,1) = TmpVal

              else
                DumpBio(i,j,k,1) = ValuesToTest(i,j,k,l)
              endif
            else
              DumpBio(i,j,k,1) = bio%InitialNut(i,j,k,l)
            endif

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


  DEALLOCATE(DumpBio, ValuesToTest)


end subroutine wrt_upd_nut
