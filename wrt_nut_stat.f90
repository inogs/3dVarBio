subroutine wrt_nut_stat

  use set_knd
  use grd_str
  use drv_str
  use mpi_str
  use bio_str
  use pnetcdf
  use da_params

  implicit none

  INTEGER(i4)        :: ncid, ierr, i, j, k, l, m, mm
  INTEGER(i4)        :: idP, iVar
  INTEGER(I4)        :: xid,yid,depid,timeId, idTim
  ! INTEGER(I4)        :: indN3n,indN1p
  INTEGER            :: system, SysErr

  INTEGER(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km, MyTime
  INTEGER(KIND=MPI_OFFSET_KIND) :: MyCountSingle(1), MyStartSingle(1)
  CHARACTER(LEN=46)    :: BioRestart
  CHARACTER(LEN=47)    :: BioRestartLong
  CHARACTER(LEN=6)     :: MyVarName

  real(r8)           :: TmpVal, MyCorr, MyRatio!, SMALL
  real(r4), allocatable, dimension(:,:,:,:) :: ValuesToTest

  ! bug fix Intel 2018
  real(r4), allocatable, dimension(:,:,:,:) :: DumpBio
!  real(r4), allocatable, dimension(:,:,:) :: Nut
  integer(KIND=MPI_OFFSET_KIND) :: MyStart_4d(4), MyCount_4d(4)

  real(r8) :: TimeArr(1)

!  SMALL = 1.e-5

  MyStart_4d(1:3) = MyStart(:)
  MyStart_4d(4) = 1
  MyCount_4d(1:3) = MyCount(:)
  MyCount_4d(4) = 1

  
  ALLOCATE(DumpBio(grd%im,grd%jm,grd%km,1)); DumpBio(:,:,:,:) = 1.e20
  ALLOCATE(ValuesToTest(grd%im,grd%jm,grd%km,NNVar)); ValuesToTest(:,:,:,:) = dble(0.)
!  ALLOCATE(Nut(grd%im,grd%jm,grd%km)); Nut(:,:,:) = dble(0.)

  if(MyId .eq. 0) then
    write(drv%dia,*) 'writing nut structure'     
  endif

  
  global_im = GlobalRow
  global_jm = GlobalCol
  global_km = grd%km
  MyTime = 1
  
  MyCountSingle(1) = 1
  MyStartSingle(1) = 1
  TimeArr(1) = DA_JulianDate
  
  
  if(bio%n3n .eq. 0) then
    write(*,*) "ERROR: Nitrate to be assimilated NOT set in namelist"
    write(drv%dia,*) "ERROR: Nitrate to be assimilated NOT set in namelist"
    call MPI_Barrier(Var3DCommunicator, ierr)
    call MPI_Abort(Var3DCommunicator,-1,ierr)
  endif
  
  if((bio%updateN1p .eq. 1) .and. (NNVar.lt.2)) then
    write(*,*) "ERROR: Required phosphate update but NOT set in DA_params.f90"
    write(drv%dia,*) "ERROR: Required phosphate update but NOT set in DA_params.f90"
    call MPI_Barrier(Var3DCommunicator, ierr)
    call MPI_Abort(Var3DCommunicator,-1,ierr)
  endif
  
  ! if(MyId .eq. 0) then
  !   write(drv%dia,*) "before 3d do"
  ! endif
  
  ! do l=1,NNutVar
  !   iVar = NPhytoVar + l
  !   if (DA_VarList(iVar) .eq. 'N3n') then
  !     indN3n = l
  !   endif
  !   if(DA_VarList(iVar) .eq. 'N1p') then
  !     indN1p = l
  !   endif
  ! enddo

  ! The following cycle is hard coded on numbers of variables --> to be updated in a more general way
  do k=1,grd%km
    do j=1,grd%jm
      do i=1,grd%im
        if(bio%InitialNut(i,j,k,1) .lt. 1.e20) then
          ! check obtained values and eventually
          ! correct them in order to avoid negative concentrations
          ! if the correction is negative, the correction must be reduced
          ! if(bio%n3n.eq.1) then
              ValuesToTest(i,j,k,1) = bio%InitialNut(i,j,k,1) + grd%n3n(i,j,k)
              if(bio%updateN1p.eq.1) then
                ValuesToTest(i,j,k,2) = bio%InitialNut(i,j,k,2) + grd%n3n(i,j,k)*bio%covn3n_n1p(i,j,k)
            !  endif
          ! else
          !    ValuesToTest(i,j,k,1) = bio%InitialNut(i,j,k,1)
          !    ValuesToTest(i,j,k,2) = bio%InitialNut(i,j,k,2)
          ! endif
          ! if(bio%o2o.eq.1) then
          !   ValuesToTest(i,j,k,3) = bio%InitialNut(i,j,k,3) + grd%o2o(i,j,k)
              endif
        endif
      enddo
    enddo
  enddo
  
  
  
  do l=1,NNVar
    iVar = NPhytoVar + l
    
    if(iVar .gt. NBioVar) then
      if(MyId .eq. 0) &
      write(*,*) "Warning: Reading a variable not in the DA_VarList!"
    endif
    
    BioRestart = 'RST_after.'//ShortDate//'.'//DA_VarList(iVar)//'.nc'
    BioRestartLong = 'RST_after.'//DA_DATE//'.'//DA_VarList(iVar)//'.nc'
    
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

            if(grd%msk(i,j,k).eq.1) then

              if(ValuesToTest(i,j,k,l) .lt. 0) then
                ! Excluding negative concentrations
                ! This correction must be the first
                ! condition applied (before apply corrections
                ! on the other components)
                TmpVal = 0.1*bio%InitialNut(i,j,k,l)
                ! if(TmpVal.gt.SMALL) then
                !   TmpVal = SMALL
                ! endif
                DumpBio(i,j,k,1) = TmpVal

              else
                DumpBio(i,j,k,1) = ValuesToTest(i,j,k,l)
                ! if(bio%ApplyConditions) then

                ! endif ! ApplyConditions

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
  ! DEALLOCATE(DumpBio, ValuesToTest, Nut)

end subroutine wrt_nut_stat
